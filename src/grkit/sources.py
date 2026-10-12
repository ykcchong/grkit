"""Uniform access to GeneReviews chapters from a directory *or* straight out of a tarball.

The NCBI archive is a single ~620 MB ``gene_NBK1116.tar.gz`` containing ~1000 ``.nxml``
chapters plus images and PDFs. Unpacking it into a directory of a thousand files (and
~870 MB on disk) is unnecessary for most work, so every consumer in this package talks to
a :class:`ChapterSource` instead:

``DirectorySource``
    An already-extracted ``nxmls/`` directory. Cheap random access.
``TarGzSource``
    A ``.tar.gz``/``.tgz``/``.tar`` archive read in place.

Two access patterns are offered because ``gzip`` streams are not seekable:

*Streaming* (:meth:`ChapterSource.iter_bytes`) decompresses exactly once and yields
chapters in archive order. This is the fast path for whole-corpus work and uses no extra
disk.

*Random access* (:meth:`ChapterSource.read`) needs a seekable stream. For a compressed
archive the first random read transparently spools the decompressed tar to a temporary
file (deleted on :meth:`close`), after which lookups are real seeks. Pass ``spool=False``
to instead rescan the stream on every read -- slower, but writes nothing.
"""

from __future__ import annotations

import io
import os
import shutil
import tarfile
import tempfile
import threading
from abc import ABC, abstractmethod
from collections.abc import Iterator
from typing import BinaryIO

NXML_SUFFIX = ".nxml"


class ChapterSource(ABC):
    """A random-access, iterable collection of ``.nxml`` chapter documents."""

    #: Human-readable description of where chapters come from.
    description: str = ""

    @abstractmethod
    def names(self) -> list[str]:
        """Sorted basenames of the available ``.nxml`` members."""

    @abstractmethod
    def read(self, name: str) -> bytes:
        """Return the raw bytes of member ``name`` (a basename)."""

    @abstractmethod
    def iter_bytes(self) -> Iterator[tuple[str, bytes]]:
        """Yield ``(name, bytes)`` for every member, in source order."""

    def open(self, name: str) -> BinaryIO:
        """Return a binary file object for ``name``."""
        return io.BytesIO(self.read(name))

    def __len__(self) -> int:
        return len(self.names())

    def __contains__(self, name: object) -> bool:
        return isinstance(name, str) and self._key(name) in set(self.names())

    @staticmethod
    def _key(name: str) -> str:
        """Normalize a lookup key to a basename, with or without the suffix."""
        base = os.path.basename(name)
        return base if base.endswith(NXML_SUFFIX) else base + NXML_SUFFIX

    def close(self) -> None:  # noqa: B027 - optional hook, not every source has resources
        """Release any held resources. Safe to call more than once."""

    def __enter__(self) -> ChapterSource:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


class DirectorySource(ChapterSource):
    """Chapters read from a directory of ``.nxml`` files."""

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = os.path.abspath(os.fspath(path))
        if not os.path.isdir(self.path):
            raise NotADirectoryError(f"not a directory: {self.path}")
        self._names: list[str] | None = None
        self.description = f"directory {self.path}"

    def names(self) -> list[str]:
        if self._names is None:
            self._names = sorted(
                entry.name
                for entry in os.scandir(self.path)
                if entry.is_file() and entry.name.endswith(NXML_SUFFIX)
            )
        return self._names

    def read(self, name: str) -> bytes:
        with open(os.path.join(self.path, self._key(name)), "rb") as handle:
            return handle.read()

    def iter_bytes(self) -> Iterator[tuple[str, bytes]]:
        for name in self.names():
            yield name, self.read(name)


class FileSource(ChapterSource):
    """A single ``.nxml`` file presented as a one-member source."""

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = os.path.abspath(os.fspath(path))
        if not os.path.isfile(self.path):
            raise FileNotFoundError(self.path)
        self._name = os.path.basename(self.path)
        self.description = f"file {self.path}"

    def names(self) -> list[str]:
        return [self._name]

    def read(self, name: str) -> bytes:
        key = self._key(name)
        if key != self._name:
            raise KeyError(f"{key!r} not in {self._name!r}")
        with open(self.path, "rb") as handle:
            return handle.read()

    def iter_bytes(self) -> Iterator[tuple[str, bytes]]:
        yield self._name, self.read(self._name)


class TarGzSource(ChapterSource):
    """Chapters read directly from a tar archive, compressed or not.

    Parameters
    ----------
    path:
        Path to the ``.tar.gz``/``.tgz``/``.tar`` archive.
    spool:
        When ``True`` (default) the first random-access read decompresses the archive to a
        temporary tar file so subsequent lookups are genuine seeks. When ``False``, each
        read rescans the compressed stream and nothing is written to disk.
    scan_budget:
        Number of random reads served by rescanning the compressed stream before falling
        back to spooling. This keeps a one-off ``get("abs")`` from writing a ~2 GB temp
        file, while a loop over many chapters still gets true random access after the
        budget is spent. Ignored when ``spool`` is ``False``.
    spool_dir:
        Directory for the temporary tar (defaults to the system temp dir).
    """

    def __init__(
        self,
        path: str | os.PathLike[str],
        *,
        spool: bool = True,
        scan_budget: int = 3,
        spool_dir: str | None = None,
    ) -> None:
        self.path = os.path.abspath(os.fspath(path))
        if not os.path.isfile(self.path):
            raise FileNotFoundError(self.path)
        self.spool = spool
        self.scan_budget = max(0, scan_budget)
        self.spool_dir = spool_dir
        self._names: list[str] | None = None
        self._spool_path: str | None = None
        self._tar: tarfile.TarFile | None = None
        self._scans = 0
        self._lock = threading.Lock()
        self.description = f"archive {self.path}"

    # ── members ──────────────────────────────────────────────────────────────
    @staticmethod
    def _is_chapter(member: tarfile.TarInfo) -> bool:
        return member.isfile() and member.name.endswith(NXML_SUFFIX)

    def names(self) -> list[str]:
        if self._names is None:
            found: list[str] = []
            with tarfile.open(self.path, "r:*") as tar:
                for member in tar:
                    if self._is_chapter(member):
                        found.append(os.path.basename(member.name))
            self._names = sorted(found)
        return self._names

    def iter_bytes(self) -> Iterator[tuple[str, bytes]]:
        """Single decompression pass over the archive, in archive order.

        Uses streaming mode (``r|*``) so memory stays flat regardless of archive size.
        """
        with tarfile.open(self.path, "r|*") as tar:
            for member in tar:
                if not self._is_chapter(member):
                    continue
                handle = tar.extractfile(member)
                if handle is None:
                    continue
                yield os.path.basename(member.name), handle.read()

    # ── random access ────────────────────────────────────────────────────────
    def _ensure_spool(self) -> tarfile.TarFile:
        with self._lock:
            if self._tar is not None:
                return self._tar
            handle = tempfile.NamedTemporaryFile(
                prefix="grkit-", suffix=".tar", dir=self.spool_dir, delete=False
            )
            self._spool_path = handle.name
            try:
                if self._looks_compressed():
                    import gzip

                    with gzip.open(self.path, "rb") as src, handle:
                        shutil.copyfileobj(src, handle, length=1 << 20)
                else:
                    with open(self.path, "rb") as src, handle:
                        shutil.copyfileobj(src, handle, length=1 << 20)
                self._tar = tarfile.open(self._spool_path, "r:")
            except Exception:
                handle.close()
                self._discard_spool()
                raise
            return self._tar

    def _looks_compressed(self) -> bool:
        with open(self.path, "rb") as handle:
            return handle.read(2) == b"\x1f\x8b"

    def _scan_for(self, key: str) -> bytes | None:
        for name, payload in self.iter_bytes():
            if name == key:
                return payload
        return None

    def read(self, name: str) -> bytes:
        key = self._key(name)
        if not self.spool:
            payload = self._scan_for(key)
            if payload is None:
                raise KeyError(f"{key!r} not in {self.path}")
            return payload

        # Serve the first few lookups by rescanning: one decompression pass costs less than
        # materialising a multi-gigabyte temp tar for a single chapter.
        if self._tar is None and self._scans < self.scan_budget:
            self._scans += 1
            payload = self._scan_for(key)
            if payload is not None:
                return payload

        tar = self._ensure_spool()
        try:
            member = tar.getmember(key)
        except KeyError:
            # Members may carry a directory prefix; fall back to a basename search.
            member = next(
                (m for m in tar.getmembers() if os.path.basename(m.name) == key), None
            )
            if member is None:
                raise KeyError(f"{key!r} not in {self.path}") from None
        handle = tar.extractfile(member)
        if handle is None:
            raise KeyError(f"{key!r} is not a regular file")
        return handle.read()

    def close(self) -> None:
        with self._lock:
            if self._tar is not None:
                self._tar.close()
                self._tar = None
            self._discard_spool()

    def _discard_spool(self) -> None:
        if self._spool_path and os.path.exists(self._spool_path):
            try:
                os.remove(self._spool_path)
            except OSError:  # pragma: no cover - best effort cleanup
                pass
        self._spool_path = None

    @property
    def spooled(self) -> bool:
        """True once the temporary decompressed tar exists."""
        return self._tar is not None


def is_archive(path: str | os.PathLike[str]) -> bool:
    """Heuristically decide whether ``path`` is a tar archive this module can read."""
    name = os.fspath(path).lower()
    if name.endswith((".tar.gz", ".tgz", ".tar")):
        return True
    if not os.path.isfile(path):
        return False
    try:
        with open(path, "rb") as handle:
            head = handle.read(2)
        if head == b"\x1f\x8b":
            return True
        return tarfile.is_tarfile(path)
    except OSError:  # pragma: no cover - unreadable file
        return False


def open_source(path: str | os.PathLike[str], **kwargs: object) -> ChapterSource:
    """Return the right :class:`ChapterSource` for ``path``.

    Accepts a directory of ``.nxml`` files, a tar archive (compressed or not), or a single
    ``.nxml`` file.
    """
    if os.path.isdir(path):
        return DirectorySource(path)
    if not os.path.exists(path):
        raise FileNotFoundError(os.fspath(path))
    if is_archive(path):
        return TarGzSource(path, **kwargs)  # type: ignore[arg-type]
    return FileSource(path)
