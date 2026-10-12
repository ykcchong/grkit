"""Data model for parsed GeneReviews chapters and their tables.

The centrepiece is :class:`Table`, which models a ``<table-wrap>`` **losslessly**:

* ``rowspan``/``colspan`` are expanded into a rectangular grid, so no cell is ever
  dropped or shifted. The old positional reader turned the canonical
  ``Gene | Method | yield`` table into rows of differing length, which silently moved a
  yield value into the wrong column whenever a gene spanned two rows.
* Multi-row ``<thead>`` structures are preserved, so transposed tables where *methods are
  columns* (``... Identified by Method`` over ``Sequence analysis`` /
  ``Gene-targeted deletion/duplication analysis``) keep their full header hierarchy.
* Footnotes keep their labels and text, and cells remember which markers they reference.

Grids are stored as ``list[list[Optional[Cell]]]`` in which a cell spanning several slots
is the *same object* at every slot it covers. :meth:`Cell.is_origin` distinguishes the
anchor slot from the continuation slots; renderers emit the anchor with ``rowspan``/
``colspan`` and skip the continuations, while :meth:`Table.to_dicts` fills the value
across the span so tabular consumers see a value in every row.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass, field


@dataclass
class Cell:
    """One table cell, with its span geometry and both text and HTML renderings."""

    text: str
    """Plain text with footnote markers removed -- what structured consumers should read."""

    html: str = ""
    """Escaped inline HTML (emphasis, breaks, links, footnote anchors)."""

    row: int = 0
    col: int = 0
    rowspan: int = 1
    colspan: int = 1
    is_header: bool = False
    scope: str = ""
    id: str = ""
    align: str = ""
    header_ids: tuple[str, ...] = ()
    """Column header ids from the ``headers`` attribute, giving semantic column mapping."""

    marks: tuple[str, ...] = ()
    """Footnote labels this cell references, all of which are defined in the table."""

    orphan_marks: tuple[str, ...] = ()
    """Marker-shaped superscripts whose footnote is missing from the source table.

    GeneReviews occasionally drops a footnote but leaves its marker in a cell (e.g. the
    ``brugada`` testing table defines 1-5, 7, 8 while a cell cites 6). These are excluded
    from :attr:`text` so parsed values stay clean, and recorded here so the gap is visible.
    """


    raw_text: str = ""
    """Plain text with footnote markers retained (useful for display/debugging)."""

    def is_origin(self, row: int, col: int) -> bool:
        """True when ``(row, col)`` is this cell's anchor slot (not a span continuation)."""
        return self.row == row and self.col == col

    @property
    def is_spanned(self) -> bool:
        return self.rowspan > 1 or self.colspan > 1

    def covers(self, row: int, col: int) -> bool:
        return (
            self.row <= row < self.row + self.rowspan
            and self.col <= col < self.col + self.colspan
        )


@dataclass
class Footnote:
    """A ``<fn>`` from ``<table-wrap-foot>``."""

    id: str
    label: str
    text: str
    html: str = ""


@dataclass
class Table:
    """A ``<table-wrap>`` reduced to label, caption, a faithful grid, and footnotes."""

    id: str = ""
    label: str = ""
    caption: str = ""
    caption_html: str = ""
    grid: list[list[Cell | None]] = field(default_factory=list)
    n_header_rows: int = 0
    footnotes: list[Footnote] = field(default_factory=list)
    section: str = ""
    chapter_id: str = ""
    n_cols: int = 0

    # ── structure ────────────────────────────────────────────────────────────
    @property
    def header_grid(self) -> list[list[Cell | None]]:
        return self.grid[: self.n_header_rows]

    @property
    def body_grid(self) -> list[list[Cell | None]]:
        return self.grid[self.n_header_rows :]

    @property
    def n_rows(self) -> int:
        return len(self.grid)

    @property
    def has_header(self) -> bool:
        return self.n_header_rows > 0

    def cells(self) -> Iterator[Cell]:
        """Yield each unique cell once, in reading order."""
        seen: set[int] = set()
        for row in self.grid:
            for cell in row:
                if cell is None or id(cell) in seen:
                    continue
                seen.add(id(cell))
                yield cell

    def origins(self) -> Iterator[Cell]:
        """Yield only anchor cells (span continuations excluded)."""
        for r, row in enumerate(self.grid):
            for c, cell in enumerate(row):
                if cell is not None and cell.is_origin(r, c):
                    yield cell

    def bucket(self, *, header: bool) -> list[list[Cell | None]]:
        return self.header_grid if header else self.body_grid

    def cell_at(self, row: int, col: int) -> Cell | None:
        if 0 <= row < len(self.grid) and 0 <= col < len(self.grid[row]):
            return self.grid[row][col]
        return None

    def column(self, index: int) -> list[Cell | None]:
        """All cells in column ``index`` (header rows first), ``None`` where empty."""
        return [row[index] if index < len(row) else None for row in self.grid]

    def data_column(self, index: int) -> list[Cell | None]:
        """Body rows only for column ``index``."""
        return self.column(index)[self.n_header_rows :]

    # ── header introspection ─────────────────────────────────────────────────
    def column_header_path(self, col: int) -> list[str]:
        """Distinct header labels covering ``col``, outermost first.

        A cell that spans several header rows (``rowspan``) contributes its label once,
        which is what makes transposed tables readable: ``Sequence analysis`` under
        ``Proportion of Pathogenic Variants Identified by Method`` yields a two-level path.
        """
        path: list[str] = []
        for row in self.header_grid:
            cell = row[col] if col < len(row) else None
            if cell is None:
                continue
            label = cell.text
            if not label or (path and path[-1] == label):
                continue
            path.append(label)
        return path

    def column_leaf_header(self, col: int) -> str:
        path = self.column_header_path(col)
        return path[-1] if path else ""

    def column_key(self, col: int) -> str:
        """A unique, human-readable dict key for column ``col``."""
        path = self.column_header_path(col)
        if not path:
            return f"col{col}"
        return " | ".join(path)

    def column_keys(self) -> list[str]:
        keys = [self.column_key(c) for c in range(self.n_cols)]
        if len(set(keys)) == len(keys):
            return keys
        # Disambiguate duplicate header paths while preserving readability.
        seen: dict[str, int] = {}
        out: list[str] = []
        for key in keys:
            seen[key] = seen.get(key, 0) + 1
            out.append(key if seen[key] == 1 else f"{key} ({seen[key]})")
        return out

    def find_columns(self, predicate: Callable[[str], bool]) -> list[int]:
        """Column indices whose leaf header satisfies ``predicate``."""
        return [c for c in range(self.n_cols) if predicate(self.column_leaf_header(c))]

    # ── tabular projections ──────────────────────────────────────────────────
    def to_grid(self, *, filled: bool = True) -> list[list[str]]:
        """Body rows as a rectangular list of string rows.

        With ``filled`` (the default) span continuations repeat the anchor value, so the
        result is aligned and rectangular; otherwise they are empty strings.
        """
        keys = self.column_keys()
        out: list[list[str]] = []
        for r in range(self.n_header_rows, len(self.grid)):
            row_out: list[str] = []
            for c in range(len(keys)):
                cell = self.grid[r][c] if c < len(self.grid[r]) else None
                if cell is None:
                    row_out.append("")
                elif filled or cell.is_origin(r, c):
                    row_out.append(cell.text)
                else:
                    row_out.append("")
            out.append(row_out)
        return out

    def to_dicts(self) -> list[dict[str, str]]:
        """Body rows as dicts keyed by header path, with spans filled in."""
        keys = self.column_keys()
        return [dict(zip(keys, row, strict=False)) for row in self.to_grid(filled=True)]

    def header_rows_text(self) -> list[list[str]]:
        """Header rows as aligned string rows (spans filled)."""
        keys = self.column_keys()
        out: list[list[str]] = []
        for r in range(self.n_header_rows):
            row_out: list[str] = []
            for c in range(len(keys)):
                cell = self.grid[r][c] if c < len(self.grid[r]) else None
                row_out.append(cell.text if cell is not None else "")
            out.append(row_out)
        return out

    def footnote_map(self) -> dict[str, Footnote]:
        return {fn.label: fn for fn in self.footnotes}

    def text(self) -> str:
        """Flattened text: label, caption, headers, body rows, then footnotes."""
        parts: list[str] = []
        if self.label:
            parts.append(self.label)
        if self.caption:
            parts.append(self.caption)
        for row in self.header_rows_text() + self.to_grid():
            parts.append(" | ".join(row))
        for fn in self.footnotes:
            parts.append(f"{fn.label} {fn.text}")
        return "\n".join(p for p in parts if p)

    # ── renderers (imported lazily to avoid a cycle) ──────────────────────────
    def to_html(self, **kwargs: object) -> str:
        from .render import table_to_html

        return table_to_html(self, **kwargs)  # type: ignore[arg-type]

    def to_markdown(self, **kwargs: object) -> str:
        from .render import table_to_markdown

        return table_to_markdown(self, **kwargs)  # type: ignore[arg-type]

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"<Table {self.id!r} {self.label!r} {self.n_rows}x{self.n_cols} "
            f"header_rows={self.n_header_rows}>"
        )


@dataclass
class Author:
    surname: str = ""
    given_names: str = ""
    email: str = ""
    affiliation: str = ""

    @property
    def full_name(self) -> str:
        return " ".join(p for p in (self.given_names, self.surname) if p)


@dataclass
class Reference:
    id: str = ""
    authors: list[str] = field(default_factory=list)
    title: str = ""
    source: str = ""
    year: str = ""
    volume: str = ""
    fpage: str = ""
    lpage: str = ""
    pmid: str = ""
    doi: str = ""
    pub_type: str = ""


@dataclass
class Link:
    text: str = ""
    href: str = ""
    link_type: str = ""


@dataclass
class CrossReference:
    """A normalized link to an external resource or another GeneReviews chapter."""

    db: str = ""
    value: str = ""
    text: str = ""
    url: str = ""


@dataclass
class Section:
    """A ``<sec>`` with its title, depth, and lazily-computed content accessors."""

    id: str = ""
    name: str = ""
    title: str = ""
    depth: int = 0
    element: object | None = field(default=None, repr=False, compare=False)
    chapter_id: str = ""

    @property
    def subsections(self) -> list[Section]:
        from .reader import sections_of

        return sections_of(self.element, depth=self.depth + 1, chapter_id=self.chapter_id)

    @property
    def paragraphs(self) -> list[str]:
        from .inline import element_to_text

        if self.element is None:
            return []
        return [
            element_to_text(p)
            for p in self.element.findall("p")  # type: ignore[attr-defined]
            if element_to_text(p)
        ]

    @property
    def tables(self) -> list[Table]:
        from .tables import tables_of

        return tables_of(self.element, section=self.name, chapter_id=self.chapter_id)

    @property
    def text(self) -> str:
        from .inline import element_to_text

        return element_to_text(self.element)  # type: ignore[arg-type]
