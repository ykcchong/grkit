"""Tests for :mod:`grkit.updater`.

Every test is hermetic: the network is replaced by fake ``urlopen`` callables and
tar.gz payloads are built in memory, so nothing here touches NCBI.
"""

import io
import json
import os
import tarfile
from urllib.error import HTTPError, URLError

import pytest

from grkit import updater

ARCHIVE_URL = f"{updater.LITARCH_BASE}/ca/84/gene_NBK1116.tar.gz"
FILE_LIST_CSV = (
    "File,Title,Publisher,Publication Year,Accession ID,Last Updated\n"
    "ab/12/gene_NBK0001.tar.gz,Some Other Book,Publisher X,2020, NBK0001 ,2020-01-01\n"
    "malformed\n"
    "ca/84/gene_NBK1116.tar.gz, GeneReviews , NCBI ,2003,NBK1116,2024-05-01 12:00:00\n"
)


class FakeResponse(io.BytesIO):
    """Minimal stand-in for the object returned by ``urllib.request.urlopen``."""

    def __init__(self, body: bytes, headers: dict[str, str] | None = None) -> None:
        super().__init__(body)
        self.headers = dict(headers or {})


class BrokenResponse(io.BytesIO):
    """Response whose body read fails midway through a download."""

    def read(self, size: int = -1) -> bytes:
        raise OSError("connection reset by peer")


def _tar_bytes(entries: list[tuple[str, bytes]]) -> bytes:
    """Build an in-memory tar.gz containing ``entries`` as regular files."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, payload in entries:
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    return buf.getvalue()


def _patch_resolve(monkeypatch, *, url: str, last_updated: str) -> None:
    """Replace :func:`updater.resolve_archive` with an offline stub."""

    def fake_resolve(
        accession: str = updater.GENEREVIEWS_ACCESSION,
        *,
        file_list_url: str = updater.FILE_LIST_URL,
        timeout: float = 120.0,
    ) -> updater.ArchiveInfo:
        return updater.ArchiveInfo(
            accession=accession,
            url=url,
            relative_path="ca/84/gene_NBK1116.tar.gz",
            title="GeneReviews",
            publisher="NCBI",
            last_updated=last_updated,
        )

    monkeypatch.setattr(updater, "resolve_archive", fake_resolve)


def _patch_resolve_error(monkeypatch, exc: Exception) -> None:
    def fake_resolve(*args: object, **kwargs: object) -> updater.ArchiveInfo:
        raise exc

    monkeypatch.setattr(updater, "resolve_archive", fake_resolve)


def _seed_cache(
    cache_dir: str,
    body: bytes,
    *,
    last_updated: str,
    url: str = ARCHIVE_URL,
    etag: str = '"v1"',
    accession: str = updater.GENEREVIEWS_ACCESSION,
) -> tuple[str, str]:
    """Write a cached tarball plus a well-formed sidecar; return their paths."""
    tar_path = os.path.join(cache_dir, f"{accession}.tar.gz")
    meta_path = os.path.join(cache_dir, f"{accession}.meta.json")
    with open(tar_path, "wb") as fh:
        fh.write(body)
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "last_updated": last_updated,
                "url": url,
                "etag": etag,
                "content_length": len(body),
                "sha256": updater.sha256_file(tar_path),
                "fetched_at": "2024-01-01T00:00:00Z",
            },
            fh,
        )
    return tar_path, meta_path


def _fail_urlopen(*args: object, **kwargs: object) -> FakeResponse:
    raise AssertionError("unexpected network access")


def _request_headers(request: object) -> dict[str, str]:
    """Lower-cased headers of a ``urllib.request.Request`` (empty for plain URLs)."""
    headers = getattr(request, "headers", {})
    return {str(key).lower(): str(value) for key, value in headers.items()}


# ──────────────────────────────────────────────────────────────────────────────
# resolve_archive
# ──────────────────────────────────────────────────────────────────────────────
def test_resolve_archive_finds_matching_row(monkeypatch) -> None:
    seen: list[str] = []

    def fake_urlopen(url: str, timeout: float | None = None) -> FakeResponse:
        seen.append(url)
        return FakeResponse(FILE_LIST_CSV.encode("utf-8"))

    monkeypatch.setattr(updater, "urlopen", fake_urlopen)

    info = updater.resolve_archive()

    assert seen == [updater.FILE_LIST_URL]
    assert info.accession == updater.GENEREVIEWS_ACCESSION
    assert info.url == ARCHIVE_URL
    assert info.relative_path == "ca/84/gene_NBK1116.tar.gz"
    assert info.title == "GeneReviews"
    assert info.publisher == "NCBI"
    assert info.last_updated == "2024-05-01 12:00:00"
    assert info.local_path == ""


def test_resolve_archive_custom_accession_and_file_list_url(monkeypatch) -> None:
    def fake_urlopen(url: str, timeout: float | None = None) -> FakeResponse:
        assert url == "https://example.test/file_list.csv"
        return FakeResponse(FILE_LIST_CSV.encode("utf-8"))

    monkeypatch.setattr(updater, "urlopen", fake_urlopen)

    info = updater.resolve_archive(
        "NBK0001", file_list_url="https://example.test/file_list.csv"
    )

    assert info.url == f"{updater.LITARCH_BASE}/ab/12/gene_NBK0001.tar.gz"
    assert info.last_updated == "2020-01-01"


def test_resolve_archive_missing_accession_raises(monkeypatch) -> None:
    monkeypatch.setattr(
        updater,
        "urlopen",
        lambda url, timeout=None: FakeResponse(FILE_LIST_CSV.encode("utf-8")),
    )

    with pytest.raises(LookupError, match="NBK9999"):
        updater.resolve_archive("NBK9999")


def test_resolve_archive_empty_file_list_raises(monkeypatch) -> None:
    monkeypatch.setattr(
        updater, "urlopen", lambda url, timeout=None: FakeResponse(b"")
    )

    with pytest.raises(LookupError, match="Empty file list"):
        updater.resolve_archive()


# ──────────────────────────────────────────────────────────────────────────────
# download_archive
# ──────────────────────────────────────────────────────────────────────────────
def test_download_archive_streams_body(monkeypatch, tmp_path) -> None:
    body = b"x" * (256 * 1024)

    def fake_urlopen(url: str, timeout: float | None = None) -> FakeResponse:
        assert url == ARCHIVE_URL
        return FakeResponse(body, {"Content-Length": str(len(body))})

    monkeypatch.setattr(updater, "urlopen", fake_urlopen)
    dest = tmp_path / "nested" / "archive.tar.gz"

    result = updater.download_archive(ARCHIVE_URL, str(dest), quiet=True)

    assert result == str(dest)
    assert dest.read_bytes() == body


def test_download_archive_without_content_length(monkeypatch, tmp_path) -> None:
    body = b"no-length-header"

    monkeypatch.setattr(
        updater, "urlopen", lambda url, timeout=None: FakeResponse(body)
    )
    dest = tmp_path / "archive.tar.gz"

    assert updater.download_archive(ARCHIVE_URL, str(dest), quiet=True) == str(dest)
    assert dest.read_bytes() == body


def test_download_archive_logs_progress(monkeypatch, tmp_path, capsys) -> None:
    body = b"y" * 64
    monkeypatch.setattr(updater, "_PROGRESS_INTERVAL", 16)  # pretend 16 bytes == 16 MB
    monkeypatch.setattr(
        updater,
        "urlopen",
        lambda url, timeout=None: FakeResponse(body, {"Content-Length": str(len(body))}),
    )

    updater.download_archive(ARCHIVE_URL, str(tmp_path / "archive.tar.gz"))

    err = capsys.readouterr().err
    assert f"Downloading {ARCHIVE_URL}" in err
    assert "100.0%" in err
    assert "Downloaded" in err


# ──────────────────────────────────────────────────────────────────────────────
# extract_nxmls
# ──────────────────────────────────────────────────────────────────────────────
def test_extract_nxmls_flattens_and_skips_non_nxml(tmp_path) -> None:
    tar_path = tmp_path / "corpus.tar.gz"
    tar_path.write_bytes(
        _tar_bytes(
            [
                ("gene_NBK1116/nested/alpha.nxml", b"<alpha/>"),
                ("gene_NBK1116/beta.nxml", b"<beta/>"),
                ("gene_NBK1116/README.txt", b"not an nxml"),
                ("gene_NBK1116/manifest.sha256", b"deadbeef"),
            ]
        )
    )
    dest = tmp_path / "out"

    extracted = updater.extract_nxmls(str(tar_path), str(dest), quiet=True)

    assert extracted == ["alpha.nxml", "beta.nxml"]
    assert (dest / "alpha.nxml").read_bytes() == b"<alpha/>"
    assert (dest / "beta.nxml").read_bytes() == b"<beta/>"
    assert not (dest / "README.txt").exists()
    assert not (dest / "manifest.sha256").exists()
    assert not (dest / "nested").exists()


def test_extract_nxmls_empty_archive_returns_empty_list(tmp_path) -> None:
    tar_path = tmp_path / "empty.tar.gz"
    tar_path.write_bytes(_tar_bytes([("only.txt", b"nope")]))
    dest = tmp_path / "out"

    assert updater.extract_nxmls(str(tar_path), str(dest), quiet=True) == []
    assert os.path.isdir(dest)


# ──────────────────────────────────────────────────────────────────────────────
# manifest helpers & diffing
# ──────────────────────────────────────────────────────────────────────────────
def test_sha256_file_and_build_manifest(tmp_path) -> None:
    (tmp_path / "a.nxml").write_bytes(b"<a/>")
    (tmp_path / "b.nxml").write_bytes(b"<b/>")

    manifest = updater.build_manifest(str(tmp_path), ["a.nxml", "b.nxml"])

    assert manifest["a.nxml"]["size"] == 4
    assert manifest["a.nxml"]["sha256"] == updater.sha256_file(str(tmp_path / "a.nxml"))
    assert manifest["a.nxml"]["sha256"] != manifest["b.nxml"]["sha256"]


def test_load_manifest_missing_file_returns_empty(tmp_path) -> None:
    assert updater.load_manifest(str(tmp_path / "nope.json")) == {"files": {}}


def test_save_and_load_manifest_roundtrip(tmp_path) -> None:
    path = tmp_path / "nested" / "manifest.json"
    manifest = {"files": {"a.nxml": {"sha256": "abc", "size": 4}}}

    updater.save_manifest(str(path), manifest)

    assert updater.load_manifest(str(path)) == manifest


def test_diff_manifests_reports_all_categories() -> None:
    old = {
        "files": {
            "same.nxml": {"sha256": "1", "size": 1},
            "changed.nxml": {"sha256": "2", "size": 2},
            "removed.nxml": {"sha256": "3", "size": 3},
        }
    }
    new = {
        "files": {
            "same.nxml": {"sha256": "1", "size": 1},
            "changed.nxml": {"sha256": "9", "size": 7},
            "added.nxml": {"sha256": "4", "size": 4},
        }
    }

    report = updater.diff_manifests(old, new)

    assert report.added == ["added.nxml"]
    assert report.modified == ["changed.nxml"]
    assert report.removed == ["removed.nxml"]
    assert report.unchanged == ["same.nxml"]
    assert report.has_changes


def test_diff_manifests_no_changes_has_changes_false() -> None:
    manifest = {"files": {"a.nxml": {"sha256": "1", "size": 1}}}

    report = updater.diff_manifests(manifest, manifest)

    assert report.unchanged == ["a.nxml"]
    assert not report.has_changes


# ──────────────────────────────────────────────────────────────────────────────
# fetch_archive
# ──────────────────────────────────────────────────────────────────────────────
def test_fetch_archive_downloads_into_cache(monkeypatch, tmp_path) -> None:
    body = _tar_bytes([("a.nxml", b"<a/>")])
    calls: list[object] = []

    def fake_urlopen(request: object, timeout: float | None = None) -> FakeResponse:
        calls.append(request)
        return FakeResponse(
            body, {"Content-Length": str(len(body)), "ETag": '"fresh"'}
        )

    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")
    monkeypatch.setattr(updater, "urlopen", fake_urlopen)

    info = updater.fetch_archive(str(tmp_path), quiet=True)

    assert isinstance(calls[0], str)  # no cache to revalidate -> plain URL request
    assert info.local_path == str(tmp_path / "NBK1116.tar.gz")
    assert (tmp_path / "NBK1116.tar.gz").read_bytes() == body

    meta = json.loads((tmp_path / "NBK1116.meta.json").read_text(encoding="utf-8"))
    assert set(meta) == {
        "last_updated",
        "url",
        "etag",
        "content_length",
        "sha256",
        "fetched_at",
    }
    assert meta["last_updated"] == "2024-05-01 12:00:00"
    assert meta["url"] == ARCHIVE_URL
    assert meta["etag"] == '"fresh"'
    assert meta["content_length"] == len(body)
    assert meta["sha256"] == updater.sha256_file(str(tmp_path / "NBK1116.tar.gz"))
    assert meta["fetched_at"].endswith("Z")


def test_fetch_archive_skips_download_when_cache_is_fresh(monkeypatch, tmp_path) -> None:
    body = _tar_bytes([("a.nxml", b"<a/>")])
    tar_path, _ = _seed_cache(
        str(tmp_path), body, last_updated="2024-05-01 12:00:00"
    )
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")
    monkeypatch.setattr(updater, "urlopen", _fail_urlopen)

    info = updater.fetch_archive(str(tmp_path), quiet=True)

    assert info.local_path == tar_path
    assert (tmp_path / "NBK1116.tar.gz").read_bytes() == body


def test_fetch_archive_redownloads_when_last_updated_changes(
    monkeypatch, tmp_path
) -> None:
    old_body = b"stale-archive"
    new_body = _tar_bytes([("b.nxml", b"<b/>")])
    _seed_cache(str(tmp_path), old_body, last_updated="2024-01-01 00:00:00")
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-06-06 06:06:06")
    calls: list[object] = []

    def fake_urlopen(request: object, timeout: float | None = None) -> FakeResponse:
        calls.append(request)
        return FakeResponse(
            new_body, {"Content-Length": str(len(new_body)), "ETag": '"v2"'}
        )

    monkeypatch.setattr(updater, "urlopen", fake_urlopen)

    info = updater.fetch_archive(str(tmp_path), quiet=True)

    assert len(calls) == 1
    assert _request_headers(calls[0])["if-none-match"] == '"v1"'
    assert info.local_path == str(tmp_path / "NBK1116.tar.gz")
    assert (tmp_path / "NBK1116.tar.gz").read_bytes() == new_body

    meta = json.loads((tmp_path / "NBK1116.meta.json").read_text(encoding="utf-8"))
    assert meta["last_updated"] == "2024-06-06 06:06:06"
    assert meta["etag"] == '"v2"'
    assert meta["sha256"] == updater.sha256_file(str(tmp_path / "NBK1116.tar.gz"))


def test_fetch_archive_force_redownloads_fresh_cache(monkeypatch, tmp_path) -> None:
    old_body = b"already-cached"
    new_body = b"forced-refresh"
    _seed_cache(str(tmp_path), old_body, last_updated="2024-05-01 12:00:00")
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")
    calls: list[object] = []

    def fake_urlopen(request: object, timeout: float | None = None) -> FakeResponse:
        calls.append(request)
        return FakeResponse(new_body, {"Content-Length": str(len(new_body))})

    monkeypatch.setattr(updater, "urlopen", fake_urlopen)

    updater.fetch_archive(str(tmp_path), force=True, quiet=True)

    assert len(calls) == 1
    assert isinstance(calls[0], str)  # forced refresh is unconditional
    assert (tmp_path / "NBK1116.tar.gz").read_bytes() == new_body


def test_fetch_archive_treats_http_304_as_current(monkeypatch, tmp_path) -> None:
    body = b"cached-archive-bytes"
    tar_path, _ = _seed_cache(
        str(tmp_path), body, last_updated="2024-01-01 00:00:00", etag='"v1"'
    )
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-07-07 07:07:07")
    calls: list[object] = []

    def fake_urlopen(request: object, timeout: float | None = None) -> FakeResponse:
        calls.append(request)
        raise HTTPError(ARCHIVE_URL, 304, "Not Modified", {}, io.BytesIO(b""))

    monkeypatch.setattr(updater, "urlopen", fake_urlopen)

    info = updater.fetch_archive(str(tmp_path), quiet=True)

    assert len(calls) == 1
    assert info.local_path == tar_path
    assert (tmp_path / "NBK1116.tar.gz").read_bytes() == body

    meta = json.loads((tmp_path / "NBK1116.meta.json").read_text(encoding="utf-8"))
    assert meta["last_updated"] == "2024-07-07 07:07:07"  # index advanced, bytes did not
    assert meta["etag"] == '"v1"'
    assert meta["sha256"] == updater.sha256_file(tar_path)
    assert not list(tmp_path.glob(".NBK1116-*"))


def test_fetch_archive_explicit_url_uses_last_modified(monkeypatch, tmp_path) -> None:
    body = _tar_bytes([("a.nxml", b"<a/>")])
    calls: list[object] = []

    def fake_urlopen(request: object, timeout: float | None = None) -> FakeResponse:
        calls.append(request)
        return FakeResponse(
            body,
            {
                "Content-Length": str(len(body)),
                "Last-Modified": "Wed, 01 May 2024 00:00:00 GMT",
            },
        )

    monkeypatch.setattr(updater, "urlopen", fake_urlopen)
    monkeypatch.setattr(
        updater, "resolve_archive", _fail_urlopen
    )  # explicit URL must not hit the index

    info = updater.fetch_archive(str(tmp_path), url=ARCHIVE_URL, quiet=True)

    assert calls == [ARCHIVE_URL]
    assert info.relative_path == "gene_NBK1116.tar.gz"
    assert info.last_updated == "Wed, 01 May 2024 00:00:00 GMT"

    meta = json.loads((tmp_path / "NBK1116.meta.json").read_text(encoding="utf-8"))
    assert meta["last_updated"] == "Wed, 01 May 2024 00:00:00 GMT"


def test_fetch_archive_accepts_minimal_sidecar(monkeypatch, tmp_path) -> None:
    """A sidecar recording only ``last_updated`` is enough to trust the cache."""
    (tmp_path / "NBK1116.tar.gz").write_bytes(b"cached-body")
    (tmp_path / "NBK1116.meta.json").write_text(
        json.dumps({"last_updated": "2024-05-01 12:00:00"}), encoding="utf-8"
    )
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")
    monkeypatch.setattr(updater, "urlopen", _fail_urlopen)

    info = updater.fetch_archive(str(tmp_path), quiet=True)

    assert info.local_path == str(tmp_path / "NBK1116.tar.gz")
    assert (tmp_path / "NBK1116.tar.gz").read_bytes() == b"cached-body"


def test_fetch_archive_cleans_up_after_failed_download(monkeypatch, tmp_path) -> None:
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")
    monkeypatch.setattr(
        updater, "urlopen", lambda request, timeout=None: BrokenResponse(b"")
    )

    with pytest.raises(OSError, match="connection reset"):
        updater.fetch_archive(str(tmp_path), quiet=True)

    assert os.listdir(str(tmp_path)) == []
    assert not (tmp_path / "NBK1116.tar.gz").exists()


def test_fetch_archive_redownloads_when_sidecar_is_missing(monkeypatch, tmp_path) -> None:
    (tmp_path / "NBK1116.tar.gz").write_bytes(b"orphaned-tarball")
    body = _tar_bytes([("a.nxml", b"<a/>")])
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")
    monkeypatch.setattr(
        updater,
        "urlopen",
        lambda request, timeout=None: FakeResponse(
            body, {"Content-Length": str(len(body))}
        ),
    )

    updater.fetch_archive(str(tmp_path), quiet=True)

    assert (tmp_path / "NBK1116.tar.gz").read_bytes() == body
    assert (tmp_path / "NBK1116.meta.json").is_file()


# ──────────────────────────────────────────────────────────────────────────────
# is_stale
# ──────────────────────────────────────────────────────────────────────────────
def test_is_stale_without_cache_is_true_and_offline(monkeypatch, tmp_path) -> None:
    _patch_resolve_error(monkeypatch, AssertionError("no network for empty cache"))

    assert updater.is_stale(str(tmp_path)) is True


def test_is_stale_false_when_sidecar_matches_remote(monkeypatch, tmp_path) -> None:
    _seed_cache(str(tmp_path), b"body", last_updated="2024-05-01 12:00:00")
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")

    assert updater.is_stale(str(tmp_path)) is False


def test_is_stale_true_when_remote_changed(monkeypatch, tmp_path) -> None:
    _seed_cache(str(tmp_path), b"body", last_updated="2024-01-01 00:00:00")
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-09-09 09:09:09")

    assert updater.is_stale(str(tmp_path)) is True


def test_is_stale_false_with_minimal_sidecar(monkeypatch, tmp_path) -> None:
    (tmp_path / "NBK1116.tar.gz").write_bytes(b"cached-body")
    (tmp_path / "NBK1116.meta.json").write_text(
        json.dumps({"last_updated": "2024-05-01 12:00:00"}), encoding="utf-8"
    )
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")

    assert updater.is_stale(str(tmp_path)) is False


def test_is_stale_returns_true_offline(monkeypatch, tmp_path) -> None:
    _seed_cache(str(tmp_path), b"body", last_updated="2024-05-01 12:00:00")
    _patch_resolve_error(monkeypatch, URLError("network is unreachable"))

    assert updater.is_stale(str(tmp_path)) is True


def test_is_stale_handles_corrupt_sidecar(monkeypatch, tmp_path) -> None:
    (tmp_path / "NBK1116.tar.gz").write_bytes(b"cached-body")
    (tmp_path / "NBK1116.meta.json").write_text("{not json", encoding="utf-8")
    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")

    assert updater.is_stale(str(tmp_path)) is True


def test_is_stale_propagates_file_list_url(monkeypatch, tmp_path) -> None:
    _seed_cache(str(tmp_path), b"body", last_updated="2024-05-01 12:00:00")
    seen: list[str] = []

    def fake_resolve(
        accession: str,
        *,
        file_list_url: str = updater.FILE_LIST_URL,
        timeout: float = 120.0,
    ) -> updater.ArchiveInfo:
        seen.append(file_list_url)
        return updater.ArchiveInfo(
            accession=accession,
            url=ARCHIVE_URL,
            relative_path="ca/84/gene_NBK1116.tar.gz",
            last_updated="2024-05-01 12:00:00",
        )

    monkeypatch.setattr(updater, "resolve_archive", fake_resolve)

    staled = updater.is_stale(str(tmp_path), file_list_url="https://mirror.test/f.csv")

    assert staled is False
    assert seen == ["https://mirror.test/f.csv"]


# ──────────────────────────────────────────────────────────────────────────────
# update_corpus
# ──────────────────────────────────────────────────────────────────────────────
def test_update_corpus_skip_download_diffs_existing_directory(tmp_path) -> None:
    nxml_dir = tmp_path / "nxml"
    nxml_dir.mkdir()
    (nxml_dir / "a.nxml").write_bytes(b"<a/>")
    (nxml_dir / "b.nxml").write_bytes(b"<b/>")
    manifest_path = tmp_path / "manifest.json"
    summary_path = tmp_path / "summary.json"

    updater.save_manifest(
        str(manifest_path),
        {
            "generated_at": "2024-01-01T00:00:00Z",
            "source_url": ARCHIVE_URL,
            "accession": updater.GENEREVIEWS_ACCESSION,
            "total_files": 2,
            "files": {
                "a.nxml": {
                    "sha256": updater.sha256_file(str(nxml_dir / "a.nxml")),
                    "size": 4,
                },
                "gone.nxml": {"sha256": "0" * 64, "size": 1},
            },
        },
    )

    report = updater.update_corpus(
        str(nxml_dir),
        str(manifest_path),
        skip_download=True,
        summary_path=str(summary_path),
        quiet=True,
    )

    assert report.added == ["b.nxml"]
    assert report.modified == []
    assert report.removed == ["gone.nxml"]
    assert report.unchanged == ["a.nxml"]
    assert report.has_changes

    new_manifest = updater.load_manifest(str(manifest_path))
    assert new_manifest["source_url"] == ARCHIVE_URL
    assert new_manifest["total_files"] == 2
    assert sorted(new_manifest["files"]) == ["a.nxml", "b.nxml"]

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["changes"]["added"] == ["b.nxml"]
    assert summary["changes"]["removed"] == ["gone.nxml"]


def test_update_corpus_downloads_extracts_and_prunes(monkeypatch, tmp_path) -> None:
    tar_body = _tar_bytes(
        [("gene_NBK1116/nested/c.nxml", b"<c/>"), ("gene_NBK1116/d.nxml", b"<d/>")]
    )
    nxml_dir = tmp_path / "nxml"
    nxml_dir.mkdir()
    (nxml_dir / "stale.nxml").write_bytes(b"<stale/>")
    manifest_path = tmp_path / "manifest.json"

    updater.save_manifest(
        str(manifest_path),
        {
            "source_url": "https://old.example/archive.tar.gz",
            "files": {
                "stale.nxml": {"sha256": "1" * 64, "size": 8},
                "c.nxml": {"sha256": "0" * 64, "size": 1},
            },
        },
    )

    def fake_download(
        url: str, dest_path: str, *, timeout: float = 600.0, quiet: bool = False
    ) -> str:
        assert url == ARCHIVE_URL
        with open(dest_path, "wb") as fh:
            fh.write(tar_body)
        return dest_path

    _patch_resolve(monkeypatch, url=ARCHIVE_URL, last_updated="2024-05-01 12:00:00")
    monkeypatch.setattr(updater, "download_archive", fake_download)

    report = updater.update_corpus(
        str(nxml_dir), str(manifest_path), prune_removed=True, quiet=True
    )

    assert report.added == ["d.nxml"]
    assert report.modified == ["c.nxml"]
    assert report.removed == ["stale.nxml"]
    assert sorted(os.listdir(nxml_dir)) == ["c.nxml", "d.nxml"]
    assert updater.load_manifest(str(manifest_path))["source_url"] == ARCHIVE_URL


def test_update_corpus_skip_download_requires_directory(tmp_path) -> None:
    with pytest.raises(NotADirectoryError):
        updater.update_corpus(
            str(tmp_path / "missing"), str(tmp_path / "manifest.json"), skip_download=True
        )
