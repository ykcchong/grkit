"""Self-updating access to the GeneReviews NXML archive on the NCBI FTP service.

GeneReviews is published as a single gzipped tarball by NCBI's Literature Archive
(``litarch``) service, but the tarball does **not** live at a stable URL. NCBI shards
the archive into a *content-hashed* directory tree (for example
``pub/litarch/ca/84/gene_NBK1116.tar.gz``) and the two hash segments change whenever
the archive is re-sharded server-side. Hard-coding that path guarantees a 404 the
first time NCBI reorganises the tree, so the current location is instead resolved
dynamically from the published ``file_list.csv`` index, keyed by the stable Bookshelf
accession ``NBK1116``.

Two layers of caching keep "update itself" cheap and idempotent:

* :func:`resolve_archive` fetches only the small CSV index and reports the current
  URL plus NCBI's ``Last Updated`` stamp for the accession.
* :func:`fetch_archive` caches the tarball *itself* under a caller-supplied directory
  next to a ``<accession>.meta.json`` sidecar. When the index still reports the same
  ``last_updated`` the cached tarball is returned untouched, so downstream code that
  can stream members straight out of the tar.gz never has to extract anything.
  A conditional ``If-None-Match`` request lets a server answer ``304 Not Modified``
  when the bytes are unchanged even though the index timestamp moved.

:func:`extract_nxmls` is therefore optional: it materialises the ``.nxml`` members
(flattened to basenames) when a plain directory of files is more convenient, and
:func:`update_corpus` wires download → extract → SHA-256 manifest → add/modify/remove
diffing into one call. :func:`is_stale` is an offline-safe predicate for callers that
only want to know whether a refresh is worth attempting.

The module is import-safe (no side effects), stdlib-only, and never required to run
in an environment that can reach the network: every network call is isolated behind a
single :func:`urllib.request.urlopen` indirection so it can be monkeypatched in tests.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import shutil
import sys
import tarfile
import tempfile
import time
from contextlib import suppress
from dataclasses import asdict, dataclass, field
from typing import NamedTuple
from urllib.error import HTTPError
from urllib.request import Request, urlopen

LITARCH_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/litarch"
FILE_LIST_URL = f"{LITARCH_BASE}/file_list.csv"
GENEREVIEWS_ACCESSION = "NBK1116"

_CHUNK_SIZE = 65536
_PROGRESS_INTERVAL = 16 * 1024 * 1024


# ──────────────────────────────────────────────────────────────────────────────
# Result objects
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class ArchiveInfo:
    """Resolved location/metadata of an NCBI archive.

    ``local_path`` is populated by :func:`fetch_archive` and points at the cached
    tarball on disk; it is empty when only the remote location is known.
    """

    accession: str
    url: str
    relative_path: str
    title: str = ""
    publisher: str = ""
    last_updated: str = ""
    local_path: str = ""


@dataclass
class ChangeReport:
    """Difference between a previous and a freshly built manifest."""

    added: list[str] = field(default_factory=list)
    modified: list[str] = field(default_factory=list)
    removed: list[str] = field(default_factory=list)
    unchanged: list[str] = field(default_factory=list)

    @property
    def has_changes(self) -> bool:
        return bool(self.added or self.modified or self.removed)


class _DownloadResult(NamedTuple):
    """Outcome of a single streamed download attempt."""

    not_modified: bool
    etag: str
    content_length: int
    last_modified: str


def _log(msg: str, *, quiet: bool = False) -> None:
    if not quiet:
        print(msg, file=sys.stderr)


def _header(resp: object, name: str) -> str:
    """Read a response header, tolerating missing headers or plain dicts."""
    headers = getattr(resp, "headers", None)
    if headers is None:
        return ""
    getter = getattr(headers, "get", None)
    if getter is None:
        return ""
    return str(getter(name) or "").strip()


def _int_or_zero(text: str) -> int:
    try:
        return int(text)
    except (TypeError, ValueError):
        return 0


def _basename_from_url(url: str) -> str:
    return url.split("?", 1)[0].rstrip("/").rsplit("/", 1)[-1]


# ──────────────────────────────────────────────────────────────────────────────
# Archive discovery
# ──────────────────────────────────────────────────────────────────────────────
def resolve_archive(
    accession: str = GENEREVIEWS_ACCESSION,
    *,
    file_list_url: str = FILE_LIST_URL,
    timeout: float = 120.0,
) -> ArchiveInfo:
    """Look up the current archive URL for ``accession`` in the NCBI file list.

    The published index is `File, Title, Publisher, Publication Year, Accession ID,
    Last Updated`; short or malformed rows are skipped rather than raising.

    Raises :class:`LookupError` if the accession is not present.
    """
    with urlopen(file_list_url, timeout=timeout) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    rows = list(csv.reader(io.StringIO(text)))
    if not rows:
        raise LookupError(f"Empty file list at {file_list_url}")

    for row in rows[1:]:
        if len(row) < 5:
            continue
        rel_path, title, publisher = row[0].strip(), row[1].strip(), row[2].strip()
        if row[4].strip() != accession:
            continue
        last_updated = row[5].strip() if len(row) > 5 else ""
        return ArchiveInfo(
            accession=accession,
            url=f"{LITARCH_BASE}/{rel_path}",
            relative_path=rel_path,
            title=title,
            publisher=publisher,
            last_updated=last_updated,
        )

    raise LookupError(
        f"Accession {accession} not found in {file_list_url} "
        f"({len(rows) - 1} entries scanned)"
    )


# ──────────────────────────────────────────────────────────────────────────────
# Download & extract
# ──────────────────────────────────────────────────────────────────────────────
def _stream_to_file(
    url: str,
    dest_path: str,
    *,
    timeout: float = 600.0,
    quiet: bool = False,
    etag: str = "",
) -> _DownloadResult:
    """Stream ``url`` into ``dest_path``, optionally revalidating with ``etag``.

    Returns a :class:`_DownloadResult`; ``not_modified`` is ``True`` (and
    ``dest_path`` is left untouched/absent) when the server answered ``304``.
    """
    request: str | Request = url
    if etag:
        request = Request(url, headers={"If-None-Match": etag})
        _log(f"Revalidating cached copy of {url}", quiet=quiet)
    else:
        _log(f"Downloading {url}", quiet=quiet)

    t0 = time.time()
    try:
        resp = urlopen(request, timeout=timeout)
    except HTTPError as exc:
        if exc.code == 304 and etag:
            _log("Server reports the cached archive is current (HTTP 304)", quiet=quiet)
            return _DownloadResult(True, etag, 0, _header(exc, "Last-Modified"))
        raise

    with resp:
        total = _int_or_zero(_header(resp, "Content-Length"))
        downloaded = 0
        next_mark = _PROGRESS_INTERVAL
        with open(dest_path, "wb") as fh:
            while True:
                chunk = resp.read(_CHUNK_SIZE)
                if not chunk:
                    break
                fh.write(chunk)
                downloaded += len(chunk)
                if downloaded >= next_mark:
                    next_mark += _PROGRESS_INTERVAL
                    pct = f"{downloaded / total * 100:5.1f}% " if total else ""
                    _log(f"  {pct}({downloaded // 1024 // 1024} MB)", quiet=quiet)
        _log(
            f"Downloaded {downloaded / 1024 / 1024:.1f} MB in {time.time() - t0:.1f}s",
            quiet=quiet,
        )
        return _DownloadResult(
            not_modified=False,
            etag=_header(resp, "ETag"),
            content_length=total or downloaded,
            last_modified=_header(resp, "Last-Modified"),
        )


def download_archive(
    url: str, dest_path: str, *, timeout: float = 600.0, quiet: bool = False
) -> str:
    """Stream ``url`` to ``dest_path``; returns ``dest_path``.

    Progress is reported to stderr every 16 MB unless ``quiet`` is set. The response
    body is streamed in 64 KB chunks so archives never have to fit in memory.
    """
    os.makedirs(os.path.dirname(os.path.abspath(dest_path)), exist_ok=True)
    _stream_to_file(url, dest_path, timeout=timeout, quiet=quiet)
    return dest_path


def extract_nxmls(tar_path: str, dest_dir: str, *, quiet: bool = False) -> list[str]:
    """Extract ``.nxml`` members (flattened to basenames) into ``dest_dir``.

    Members carrying directory prefixes are flattened, and anything that is not a
    regular ``.nxml`` file (directories, READMEs, checksums) is skipped. Returns the
    sorted list of extracted filenames.
    """
    _log(f"Extracting .nxml files to {dest_dir}", quiet=quiet)
    os.makedirs(dest_dir, exist_ok=True)
    filenames: list[str] = []
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar.getmembers():
            if not (member.isfile() and member.name.endswith(".nxml")):
                continue
            basename = os.path.basename(member.name)
            if not basename:
                continue
            src = tar.extractfile(member)
            if src is None:
                continue
            with src, open(os.path.join(dest_dir, basename), "wb") as out:
                shutil.copyfileobj(src, out)
            filenames.append(basename)
    _log(f"Extracted {len(filenames)} .nxml files", quiet=quiet)
    return sorted(filenames)


# ──────────────────────────────────────────────────────────────────────────────
# Manifest & diffing
# ──────────────────────────────────────────────────────────────────────────────
def sha256_file(path: str) -> str:
    """Return the hex SHA-256 digest of the file at ``path``."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()


def build_manifest(nxml_dir: str, filenames: list[str]) -> dict:
    """Hash ``filenames`` (relative to ``nxml_dir``) into a manifest ``files`` map."""
    files = {}
    for name in filenames:
        path = os.path.join(nxml_dir, name)
        files[name] = {"sha256": sha256_file(path), "size": os.path.getsize(path)}
    return files


def load_manifest(path: str) -> dict:
    """Load a manifest, returning an empty one when ``path`` does not exist."""
    if os.path.isfile(path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    return {"files": {}}


def save_manifest(path: str, manifest: dict) -> None:
    """Write ``manifest`` as pretty-printed UTF-8 JSON, creating parent dirs."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def diff_manifests(old: dict, new: dict) -> ChangeReport:
    """Compare two manifests by SHA-256 and report added/modified/removed/unchanged."""
    old_files = old.get("files", {})
    new_files = new.get("files", {})
    report = ChangeReport()
    for name, info in new_files.items():
        if name not in old_files:
            report.added.append(name)
        elif old_files[name].get("sha256") != info.get("sha256"):
            report.modified.append(name)
        else:
            report.unchanged.append(name)
    report.removed = [n for n in old_files if n not in new_files]
    for lst in (report.added, report.modified, report.removed, report.unchanged):
        lst.sort()
    return report


# ──────────────────────────────────────────────────────────────────────────────
# Self-update cache
# ──────────────────────────────────────────────────────────────────────────────
def _load_meta(path: str) -> dict:
    """Load a cache sidecar, returning ``{}`` when absent, unreadable or corrupt."""
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _save_meta(path: str, meta: dict) -> None:
    """Atomically write a cache sidecar so readers never see a half-written file."""
    tmp = f"{path}.tmp"
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, ensure_ascii=False)
            fh.write("\n")
        os.replace(tmp, path)
    finally:
        with suppress(OSError):
            os.remove(tmp)


def _write_meta(
    meta_path: str,
    info: ArchiveInfo,
    *,
    last_updated: str,
    etag: str,
    content_length: int,
    sha256: str,
) -> None:
    _save_meta(
        meta_path,
        {
            "last_updated": last_updated,
            "url": info.url,
            "etag": etag,
            "content_length": content_length,
            "sha256": sha256,
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
    )


def _cache_matches(meta: dict, info: ArchiveInfo, tar_path: str) -> bool:
    """True when the sidecar proves the cached tarball is still current."""
    if not meta or not os.path.isfile(tar_path):
        return False
    recorded_url = meta.get("url")
    if recorded_url and recorded_url != info.url:
        return False
    if info.last_updated:
        return str(meta.get("last_updated") or "") == info.last_updated
    # Freshness of the remote object is unknown (explicit URL): the sidecar for this
    # exact URL is the best available evidence, so trust it and skip the download.
    return True


def fetch_archive(
    cache_dir: str,
    *,
    accession: str = GENEREVIEWS_ACCESSION,
    url: str | None = None,
    force: bool = False,
    timeout: float = 600.0,
    quiet: bool = False,
) -> ArchiveInfo:
    """Return metadata for the archive, caching the tarball inside ``cache_dir``.

    The tarball is stored at ``{cache_dir}/{accession}.tar.gz`` with a metadata sidecar
    at ``{cache_dir}/{accession}.meta.json`` recording ``last_updated``, ``url``,
    ``etag``, ``content_length``, ``sha256`` and ``fetched_at``.

    When the cache already matches the remote ``last_updated`` (and ``force`` is not
    set) nothing is downloaded and the returned :class:`ArchiveInfo` simply carries
    ``local_path``. Otherwise the body is streamed to a temporary file in ``cache_dir``
    and atomically moved into place; a partial download is always removed on failure.
    If a cached tarball exists, the request is conditional (``If-None-Match``) and an
    HTTP ``304`` marks the cached bytes as current.

    Parameters
    ----------
    url:
        Explicit archive URL. When ``None`` the URL is resolved from the NCBI index
        via ``accession``.
    force:
        Ignore the cache and unconditionally re-download.
    """
    os.makedirs(cache_dir, exist_ok=True)

    if url:
        info = ArchiveInfo(
            accession=accession,
            url=url,
            relative_path=_basename_from_url(url),
        )
    else:
        info = resolve_archive(accession)

    tar_path = os.path.join(cache_dir, f"{accession}.tar.gz")
    meta_path = os.path.join(cache_dir, f"{accession}.meta.json")
    meta = _load_meta(meta_path)

    if not force and _cache_matches(meta, info, tar_path):
        _log(
            f"Using cached archive {tar_path} (updated {info.last_updated or 'unknown'})",
            quiet=quiet,
        )
        info.local_path = tar_path
        return info

    # Only revalidate against the ETag of the same URL we recorded; a re-sharded path
    # is a different resource and must be fetched unconditionally.
    etag = ""
    if not force and os.path.isfile(tar_path) and meta.get("url") == info.url:
        etag = str(meta.get("etag") or "")

    tmp = tempfile.NamedTemporaryFile(
        dir=cache_dir, prefix=f".{accession}-", suffix=".tar.gz.part", delete=False
    )
    tmp_path = tmp.name
    tmp.close()
    try:
        result = _stream_to_file(
            info.url, tmp_path, timeout=timeout, quiet=quiet, etag=etag
        )
        if result.not_modified:
            _write_meta(
                meta_path,
                info,
                last_updated=info.last_updated or str(meta.get("last_updated") or ""),
                etag=result.etag or etag,
                content_length=_int_or_zero(str(meta.get("content_length") or "")),
                sha256=str(meta.get("sha256") or sha256_file(tar_path)),
            )
            info.local_path = tar_path
            return info

        digest = sha256_file(tmp_path)
        os.replace(tmp_path, tar_path)
        last_updated = info.last_updated or result.last_modified
        _write_meta(
            meta_path,
            info,
            last_updated=last_updated,
            etag=result.etag,
            content_length=result.content_length,
            sha256=digest,
        )
        if not info.last_updated:
            info.last_updated = last_updated
        info.local_path = tar_path
        _log(f"Cached archive at {tar_path}", quiet=quiet)
        return info
    finally:
        # Covers the failure paths and the 304 case; after os.replace the temp path is
        # already gone, and remove() then raises FileNotFoundError which we suppress.
        with suppress(OSError):
            os.remove(tmp_path)


def is_stale(
    cache_dir: str,
    *,
    accession: str = GENEREVIEWS_ACCESSION,
    file_list_url: str = FILE_LIST_URL,
) -> bool:
    """Report whether the cached tarball should be re-fetched.

    Returns ``True`` when the cache or its sidecar is missing, or when the remote
    ``last_updated`` differs from the recorded one. Any network or OS failure yields
    ``True`` as well (fail-open toward refreshing) and is never raised, so the check is
    safe to call offline or from a hot path.
    """
    tar_path = os.path.join(cache_dir, f"{accession}.tar.gz")
    meta_path = os.path.join(cache_dir, f"{accession}.meta.json")
    if not (os.path.isfile(tar_path) and os.path.isfile(meta_path)):
        return True

    meta = _load_meta(meta_path)

    try:
        info = resolve_archive(accession, file_list_url=file_list_url)
    except Exception:  # offline checks must never raise; fail open toward refresh
        return True

    return str(meta.get("last_updated") or "") != info.last_updated


# ──────────────────────────────────────────────────────────────────────────────
# Orchestrator
# ──────────────────────────────────────────────────────────────────────────────
def update_corpus(
    nxml_dir: str,
    manifest_path: str,
    *,
    url: str | None = None,
    accession: str = GENEREVIEWS_ACCESSION,
    skip_download: bool = False,
    prune_removed: bool = False,
    summary_path: str | None = None,
    quiet: bool = False,
) -> ChangeReport:
    """Refresh ``nxml_dir`` from NCBI and update ``manifest_path``.

    Parameters
    ----------
    url:
        Explicit archive URL. When ``None`` the URL is resolved dynamically from the
        NCBI file list using ``accession``.
    skip_download:
        Re-build the manifest/diff from the existing ``nxml_dir`` without downloading.
    prune_removed:
        Delete local ``.nxml`` files that are no longer present in the new archive.
    summary_path:
        Optional path for a machine-readable build summary JSON.
    """
    nxml_dir = os.path.abspath(nxml_dir)
    manifest_path = os.path.abspath(manifest_path)
    old_manifest = load_manifest(manifest_path)
    source_url = url or ""

    if skip_download:
        _log("Skipping download — using existing corpus", quiet=quiet)
        if not os.path.isdir(nxml_dir):
            raise NotADirectoryError(f"NXML directory not found: {nxml_dir}")
        filenames = sorted(f for f in os.listdir(nxml_dir) if f.endswith(".nxml"))
        source_url = old_manifest.get("source_url", "")
    else:
        if not url:
            info = resolve_archive(accession)
            source_url = info.url
            _log(
                f"Resolved {accession} -> {info.relative_path} "
                f"(updated {info.last_updated})",
                quiet=quiet,
            )
        else:
            source_url = url

        with tempfile.TemporaryDirectory() as tmp:
            tar_path = os.path.join(tmp, "genereviews.tar.gz")
            download_archive(source_url, tar_path, quiet=quiet)
            staging = os.path.join(tmp, "extracted")
            filenames = extract_nxmls(tar_path, staging, quiet=quiet)

            os.makedirs(nxml_dir, exist_ok=True)
            for name in filenames:
                shutil.move(
                    os.path.join(staging, name), os.path.join(nxml_dir, name)
                )

    new_manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source_url": source_url,
        "accession": accession,
        "total_files": len(filenames),
        "files": build_manifest(nxml_dir, filenames),
    }

    report = diff_manifests(old_manifest, new_manifest)

    if prune_removed and report.removed and not skip_download:
        for name in report.removed:
            target = os.path.join(nxml_dir, name)
            if os.path.isfile(target):
                os.remove(target)
        _log(f"Pruned {len(report.removed)} removed file(s)", quiet=quiet)

    save_manifest(manifest_path, new_manifest)
    _log(f"Manifest saved to {manifest_path}", quiet=quiet)

    if summary_path:
        summary = {
            "timestamp": new_manifest["generated_at"],
            "source_url": source_url,
            "accession": accession,
            "total_files": new_manifest["total_files"],
            "changes": asdict(report),
        }
        os.makedirs(os.path.dirname(os.path.abspath(summary_path)), exist_ok=True)
        with open(summary_path, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=2, ensure_ascii=False)
        _log(f"Build summary saved to {summary_path}", quiet=quiet)

    return report
