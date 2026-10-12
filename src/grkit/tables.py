"""Parse ``<table-wrap>`` elements into the span-aware :class:`~grkit.model.Table`.

GeneReviews tables are heavily spanned: across the corpus, 697 of the 868 canonical
"testing used in ..." tables use ``rowspan`` in the body and 178 use spans in the header,
and 172 have a two- or three-row ``<thead>``. Any parser that reads cells positionally
per ``<tr>`` will misalign those tables -- most damagingly the
``Gene | Method | Proportion ... Detectable by Method`` table, where a gene spanning two
methods makes the second row shorter and shifts the yield into the ``Method`` column.

This module expands spans into a rectangular grid instead, and keeps header hierarchy,
footnote labels, and inline markup intact.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

from .inline import (
    InlineOptions,
    element_to_html,
    element_to_text,
    element_to_text_and_marks,
    localname,
)
from .model import Cell, Footnote, Table

_CELL_TAGS = ("td", "th")


def _int_attr(el: ET.Element, name: str, default: int = 1) -> int:
    """Read a positive integer attribute, tolerating junk and clamping to ``>= 1``."""
    raw = (el.get(name) or "").strip()
    if not raw:
        return default
    try:
        value = int(raw)
    except ValueError:
        return default
    return value if value >= 1 else 1


def _cells_of(tr: ET.Element) -> list[ET.Element]:
    return [child for child in tr if localname(child.tag) in _CELL_TAGS]


def _normalize_footnote_label(raw: str) -> str:
    return raw.strip().rstrip(".):").strip()


def _footnotes(tw: ET.Element, opts_anchor: str) -> list[Footnote]:
    """Collect ``<fn>`` definitions from ``<table-wrap-foot>``.

    GeneReviews mixes numbered footnotes with *unlabelled* ones that act as a legend
    (``NA = not applicable``). Unlabelled footnotes keep an empty label rather than being
    numbered by position, which would collide with a real ``1.`` label and produce
    duplicate anchor ids.
    """
    foot = tw.find("table-wrap-foot")
    fn_els: list[ET.Element] = []
    if foot is not None:
        fn_els = [fn for fn in foot.findall("fn")]
        if not fn_els:
            fn_els = [fn for fn in foot.findall(".//fn")]
    if not fn_els:
        return []

    out: list[Footnote] = []
    for fn in fn_els:
        label_el = fn.find("label")
        label = _normalize_footnote_label(element_to_text(label_el)) if label_el is not None else ""
        # Prefer paragraph content so the footnote's own label is not repeated in the text.
        paragraphs = fn.findall("p")
        if paragraphs:
            text = " ".join(t for t in (element_to_text(p) for p in paragraphs) if t)
            html = "".join(
                element_to_html(p, InlineOptions(anchor_prefix=opts_anchor)) for p in paragraphs
            )
        else:
            text = element_to_text(fn)
            html = element_to_html(fn, InlineOptions(anchor_prefix=opts_anchor))
        out.append(Footnote(id=fn.get("id", ""), label=label, text=text.strip(), html=html))
    return out


def _split_rows(table: ET.Element) -> tuple[list[ET.Element], list[ET.Element]]:
    """Return ``(header_rows, body_rows)``.

    Rows come from ``<thead>``/``<tbody>`` when present. When a table has no ``<thead>``
    but opens with an all-``<th>`` row, that row is promoted to a header row -- unless its
    cells are ``scope="row"``, which marks body row labels rather than column headings.
    """
    thead = table.find("thead")
    tbody = table.find("tbody")

    header_els: list[ET.Element] = list(thead.findall("tr")) if thead is not None else []
    if tbody is not None:
        body_els = list(tbody.findall("tr"))
    else:
        body_els = [tr for tr in table.findall("tr")]

    if not header_els:
        promote = 0
        for tr in body_els:
            cells = _cells_of(tr)
            if not cells or not all(localname(c.tag) == "th" for c in cells):
                break
            if any((c.get("scope") or "") == "row" for c in cells):
                break
            promote += 1
        header_els = body_els[:promote]
        body_els = body_els[promote:]

    return header_els, body_els


def _make_cell(
    el: ET.Element,
    row: int,
    col: int,
    *,
    is_header: bool,
    opts: InlineOptions,
) -> Cell:
    orphans: list[str] = []
    text, marks = element_to_text_and_marks(el, opts.footnote_labels, orphans)
    return Cell(
        text=text,
        html=element_to_html(el, opts),
        row=row,
        col=col,
        rowspan=_int_attr(el, "rowspan"),
        colspan=_int_attr(el, "colspan"),
        is_header=is_header,
        scope=el.get("scope", ""),
        id=el.get("id", ""),
        align=el.get("align", ""),
        header_ids=tuple((el.get("headers") or "").split()),
        marks=marks,
        orphan_marks=tuple(orphans),
        raw_text=element_to_text(el),
    )


def _build_grid(
    rows: list[ET.Element],
    header_count: int,
    opts: InlineOptions,
) -> list[list[Cell | None]]:
    """Expand ``rows`` into a rectangular grid, honouring ``rowspan``/``colspan``.

    A spanning cell is placed as the *same object* at every slot it covers; the anchor is
    at ``(cell.row, cell.col)``. Pending downward spans are consumed before each new cell
    so columns never drift.
    """
    grid: list[list[Cell | None]] = []
    pending: dict[int, list] = {}  # column -> [rows_remaining, Cell]

    for r, tr in enumerate(rows):
        source = _cells_of(tr)
        row: list[Cell | None] = []
        x = 0
        index = 0
        while True:
            if x in pending:
                remaining, cell = pending[x]
                row.append(cell)
                if remaining <= 1:
                    del pending[x]
                else:
                    pending[x] = [remaining - 1, cell]
                x += 1
                continue
            if index < len(source):
                el = source[index]
                index += 1
                is_header_row = r < header_count
                cell = _make_cell(
                    el,
                    r,
                    x,
                    is_header=is_header_row or localname(el.tag) == "th",
                    opts=opts,
                )
                # A header cell groups *header* rows; it must never extend into the body.
                # GeneReviews emits stray `rowspan="2"` on single-row <thead> cells, which
                # would otherwise swallow the first data row and shift every column.
                if is_header_row:
                    cell.rowspan = min(cell.rowspan, header_count - r)
                for dx in range(cell.colspan):
                    row.append(cell)
                    if cell.rowspan > 1:
                        pending[x + dx] = [cell.rowspan - 1, cell]
                x += cell.colspan
                continue
            break
        grid.append(row)

    width = max((len(r) for r in grid), default=0)
    for row in grid:
        row.extend([None] * (width - len(row)))
    return grid


def parse_table_wrap(
    tw: ET.Element, *, chapter_id: str = "", section: str = ""
) -> Table:
    """Parse a single ``<table-wrap>`` element into a :class:`Table`."""
    table_id = tw.get("id", "")
    label = element_to_text(tw.find("label"))
    caption_el = tw.find("caption")
    caption = element_to_text(caption_el)

    anchor = f"fn-{table_id}" if table_id else "fn"
    footnotes = _footnotes(tw, anchor)
    # Unlabelled legend footnotes (empty label) cannot be referenced by a marker.
    labels = frozenset(fn.label for fn in footnotes if fn.label)
    opts = InlineOptions(footnote_labels=labels, anchor_prefix=anchor)
    caption_html = element_to_html(caption_el, opts)

    table = tw.find("table")
    if table is None:
        return Table(
            id=table_id,
            label=label,
            caption=caption,
            caption_html=caption_html,
            footnotes=footnotes,
            section=section,
            chapter_id=chapter_id,
            n_cols=0,
        )

    header_els, body_els = _split_rows(table)
    grid = _build_grid(header_els + body_els, len(header_els), opts)
    n_cols = max((len(row) for row in grid), default=0)

    return Table(
        id=table_id,
        label=label,
        caption=caption,
        caption_html=caption_html,
        grid=grid,
        n_header_rows=len(header_els),
        footnotes=footnotes,
        section=section,
        chapter_id=chapter_id,
        n_cols=n_cols,
    )


def tables_of(
    scope: ET.Element | None, *, chapter_id: str = "", section: str = ""
) -> list[Table]:
    """Parse every ``<table-wrap>`` nested inside ``scope``."""
    if scope is None:
        return []
    return [
        parse_table_wrap(tw, chapter_id=chapter_id, section=section)
        for tw in scope.findall(".//table-wrap")
    ]


def table_to_record(table: Table, *, include_html: bool = True) -> dict:
    """Serialize ``table`` to a JSON-friendly dict, preserving layout.

    Keys: ``id``, ``label``, ``caption``, ``section``, ``chapter_id``, ``n_rows``,
    ``n_cols``, ``n_header_rows``, ``columns`` (composite header keys), ``header_rows``,
    ``rows`` (dicts), ``grid`` (rectangular, spans filled), ``footnotes``, and optionally
    ``html`` -- the faithful rendering with real ``rowspan``/``colspan``.
    """
    record: dict = {
        "id": table.id,
        "label": table.label,
        "caption": table.caption,
        "section": table.section,
        "chapter_id": table.chapter_id,
        "n_rows": table.n_rows,
        "n_cols": table.n_cols,
        "n_header_rows": table.n_header_rows,
        "columns": table.column_keys(),
        "header_rows": table.header_rows_text(),
        "rows": table.to_dicts(),
        "grid": table.to_grid(filled=True),
        "footnotes": [
            {"id": fn.id, "label": fn.label, "text": fn.text} for fn in table.footnotes
        ],
    }
    if include_html:
        record["html"] = table.to_html()
    return record
