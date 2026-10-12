"""Corpus-level access: iterate, look up, search, and extract across all chapters.

:class:`Corpus` wraps a :class:`~grkit.sources.ChapterSource` (a directory, a ``.tar.gz``,
or a single file) and adds lazy parsing, an optional chapter cache, gene/disease lookup,
and shortcuts to the table and testing-strategy extractors.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Iterator

from .reader import Chapter
from .sources import NXML_SUFFIX, ChapterSource, open_source


class Corpus:
    """A collection of GeneReviews chapters backed by a chapter source."""

    def __init__(
        self,
        source: ChapterSource | str | os.PathLike[str],
        *,
        cache: bool = True,
        max_cache: int = 256,
    ) -> None:
        if isinstance(source, ChapterSource):
            self.source = source
            self._owns_source = False
        else:
            self.source = open_source(source)
            self._owns_source = True
        self._cache: dict[str, Chapter] = {}
        self._use_cache = cache
        self._max_cache = max_cache
        self._by_id: dict[str, str] | None = None

    # ── membership ───────────────────────────────────────────────────────────
    @property
    def names(self) -> list[str]:
        """Sorted ``.nxml`` member names."""
        return self.source.names()

    @property
    def ids(self) -> list[str]:
        """Sorted chapter ids (filename stems)."""
        return [n[: -len(NXML_SUFFIX)] for n in self.names]

    def __len__(self) -> int:
        return len(self.names)

    def __contains__(self, chapter_id: object) -> bool:
        return isinstance(chapter_id, str) and chapter_id in self._id_index()

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"<Corpus {len(self)} chapters from {self.source.description}>"

    def _id_index(self) -> dict[str, str]:
        if self._by_id is None:
            self._by_id = {n[: -len(NXML_SUFFIX)]: n for n in self.names}
        return self._by_id

    def _chapter_name(self, chapter_id: str) -> str | None:
        index = self._id_index()
        if chapter_id in index:
            return index[chapter_id]
        candidate = chapter_id if chapter_id.endswith(NXML_SUFFIX) else chapter_id + NXML_SUFFIX
        return candidate if candidate in set(self.names) else index.get(
            os.path.basename(chapter_id).removesuffix(NXML_SUFFIX)
        )

    # ── access ───────────────────────────────────────────────────────────────
    def get(self, chapter_id: str, *, raise_on_error: bool = False) -> Chapter | None:
        """Return the chapter with ``chapter_id``, or ``None`` if absent/unparseable."""
        key = os.path.basename(chapter_id).removesuffix(NXML_SUFFIX)
        if self._use_cache and key in self._cache:
            return self._cache[key]
        name = self._chapter_name(chapter_id)
        if name is None:
            return None
        try:
            chapter = Chapter.from_file(self.source.open(name))
        except Exception:
            if raise_on_error:
                raise
            return None
        chapter.source_path = name  # type: ignore[misc]
        if self._use_cache:
            if len(self._cache) >= self._max_cache:
                self._cache.pop(next(iter(self._cache)))
            self._cache[key] = chapter
        return chapter

    def iter_chapters(
        self, *, skip_errors: bool = True, names: list[str] | None = None
    ) -> Iterator[Chapter]:
        """Parse and yield every chapter, streaming from the source.

        ``skip_errors`` (the default) silently skips members with malformed XML; set it to
        ``False`` to surface parse failures.
        """
        wanted = set(names) if names is not None else None
        for name, payload in self.source.iter_bytes():
            if wanted is not None and name not in wanted:
                continue
            try:
                chapter = Chapter.from_string(payload)
            except Exception:
                if skip_errors:
                    continue
                raise
            chapter.source_path = name  # type: ignore[misc]
            yield chapter

    def __iter__(self) -> Iterator[Chapter]:
        return self.iter_chapters()

    # ── selection ────────────────────────────────────────────────────────────
    def find_by_gene(self, symbol: str, *, limit: int | None = None) -> list[Chapter]:
        """Chapters whose indexed gene symbols include ``symbol`` (case-insensitive)."""
        return self._find(lambda ch: symbol.lower() in {g.lower() for g in ch.genes}, limit)

    def find_by_disease(self, name: str, *, limit: int | None = None) -> list[Chapter]:
        """Chapters matching a disease name, in indexed or title text."""
        needle = name.lower()
        return self._find(
            lambda ch: needle in ch.title.lower()
            or any(needle in d.lower() for d in ch.diseases + ch.included_phenotypes),
            limit,
        )

    def search(self, needle: str, *, limit: int | None = None) -> list[Chapter]:
        """Chapters whose title, genes, diseases, synonyms, or body text contain ``needle``."""
        low = needle.lower()
        return self._find(
            lambda ch: low in ch.title.lower()
            or any(low in g.lower() for g in ch.genes)
            or any(low in d.lower() for d in ch.diseases + ch.synonyms)
            or low in ch.body_text.lower(),
            limit,
        )

    def _find(self, predicate: Callable[[Chapter], bool], limit: int | None) -> list[Chapter]:
        out: list[Chapter] = []
        for chapter in self.iter_chapters():
            if predicate(chapter):
                out.append(chapter)
                if limit is not None and len(out) >= limit:
                    break
        return out

    # ── tables ───────────────────────────────────────────────────────────────
    def iter_tables(self) -> Iterator[tuple[Chapter, object]]:
        """Yield ``(chapter, table)`` for every table in the corpus."""
        for chapter in self.iter_chapters():
            for table in chapter.tables:
                yield chapter, table

    def iter_testing_yield(self) -> Iterator[object]:
        """Yield a :class:`~grkit.testing.TestingYieldTable` for every testing table."""
        from .testing import extract_testing_yield

        for chapter in self.iter_chapters():
            yield from extract_testing_yield(chapter)

    def close(self) -> None:
        if self._owns_source:
            self.source.close()

    def __enter__(self) -> Corpus:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
