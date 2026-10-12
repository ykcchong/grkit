"""Render :class:`~grkit.model.Table` objects to HTML and Markdown.

The HTML output is *layout-faithful*: a cell that spans two rows is emitted once with
``rowspan="2"`` and its continuation slots are skipped, so the rendered result matches the
published GeneReviews table exactly. Header cells keep their ``scope``, multi-row headers
keep their hierarchy, and footnote markers in cells link to the footnote list rendered in
``<tfoot>``.

Cell and caption content is treated as pre-escaped HTML produced by
:mod:`grkit.inline`; everything else (ids, classes, hrefs) is escaped here.
"""

from __future__ import annotations

from collections.abc import Iterable

from .model import Footnote, Table


def footnote_anchor(prefix: str, label: str) -> str:
    """Anchor id for a footnote, matching what :mod:`grkit.inline` links to."""
    return f"{prefix}-{label}" if prefix else label


def _anchor_prefix(table: Table, id_prefix: str) -> str:
    base = f"{id_prefix}{table.id}" if table.id else ""
    return f"fn-{base}" if base else "fn"


def _attr_text(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _rewrite_anchors(html: str, old_prefix: str, new_prefix: str) -> str:
    """Re-point footnote links at a prefixed anchor id.

    Cell HTML is produced at parse time, when only the table's own id is known. When a
    caller renders several chapters into one document with ``id_prefix``, the marker links
    must follow the prefixed footnote ids or they dangle.
    """
    if not html or not old_prefix or old_prefix == new_prefix:
        return html
    return html.replace(f'href="#{old_prefix}-', f'href="#{new_prefix}-')


def _cell_attrs(cell, id_prefix: str) -> str:
    attrs: list[str] = []
    if cell.id:
        attrs.append(f'id="{_attr_text(id_prefix + cell.id)}"')
    if cell.rowspan > 1:
        attrs.append(f'rowspan="{cell.rowspan}"')
    if cell.colspan > 1:
        attrs.append(f'colspan="{cell.colspan}"')
    if cell.is_header and cell.scope:
        attrs.append(f'scope="{_attr_text(cell.scope)}"')
    if cell.header_ids:
        # Header ids carry the same prefix, so the references must too.
        attrs.append(
            f'headers="{" ".join(_attr_text(id_prefix + h) for h in cell.header_ids)}"'
        )
    if cell.align:
        attrs.append(f'style="text-align:{_attr_text(cell.align)}"')
    return (" " + " ".join(attrs)) if attrs else ""


def _render_rows(
    table: Table,
    row_range: range,
    *,
    id_prefix: str,
    old_prefix: str = "",
    new_prefix: str = "",
) -> str:
    lines: list[str] = []
    for r in row_range:
        cells = table.grid[r] if r < len(table.grid) else []
        lines.append("  <tr>")
        for c, cell in enumerate(cells):
            if cell is None or not cell.is_origin(r, c):
                continue  # continuation slot of a spanning cell: already emitted
            name = "th" if cell.is_header else "td"
            content = _rewrite_anchors(cell.html or _attr_text(cell.text), old_prefix, new_prefix)
            lines.append(f"    <{name}{_cell_attrs(cell, id_prefix)}>{content}</{name}>")
        lines.append("  </tr>")
    return "\n".join(lines)


def _render_footnotes(table: Table, prefix: str) -> str:
    if not table.footnotes:
        return ""
    lines = [
        "  <tfoot>",
        f'    <tr><td colspan="{max(table.n_cols, 1)}">',
        '      <ol class="gr-footnotes">',
    ]
    for fn in table.footnotes:
        body = fn.html or _attr_text(fn.text)
        if not fn.label:
            # Unlabelled legend note, e.g. "NA = not applicable": not a numbered footnote.
            lines.append(f'        <li class="gr-fn-legend">{body}</li>')
            continue
        anchor = footnote_anchor(prefix, fn.label)
        lines.append(
            f'        <li id="{_attr_text(anchor)}" value="{_attr_text(fn.label)}">'
            f'<span class="fn-label">{_attr_text(fn.label)}</span> {body}</li>'
        )
    lines.extend(["      </ol>", "    </td></tr>", "  </tfoot>"])
    return "\n".join(lines)


def table_to_html(
    table: Table,
    *,
    include_caption: bool = True,
    include_footnotes: bool = True,
    css_class: str = "gr-table",
    id_prefix: str = "",
    caption_class: str = "gr-caption",
) -> str:
    """Render ``table`` as an HTML ``<table>`` preserving its original layout.

    Parameters
    ----------
    include_caption:
        Emit ``<caption>`` with the table label and caption text.
    include_footnotes:
        Emit ``<tfoot>`` with the footnote definitions, targeted by the in-cell anchors.
    id_prefix:
        Prepended to the table's own id and to cell/footnote ids, so several chapters can
        be embedded in one HTML document without anchor collisions.
    """
    prefix = _anchor_prefix(table, id_prefix)
    old_prefix = f"fn-{table.id}" if table.id else "fn"
    attrs = f' class="{_attr_text(css_class)}"'
    if table.id:
        attrs += f' id="{_attr_text(id_prefix + table.id)}"'

    parts: list[str] = [f"<table{attrs}>"]

    if include_caption and (table.label or table.caption):
        caption = table.caption_html or _attr_text(table.caption)
        label = f'<span class="gr-label">{_attr_text(table.label)}</span> ' if table.label else ""
        parts.append(f'  <caption class="{_attr_text(caption_class)}">{label}{caption}</caption>')

    if table.n_header_rows:
        parts.append("  <thead>")
        parts.append(
            _render_rows(
                table,
                range(0, table.n_header_rows),
                id_prefix=id_prefix,
                old_prefix=old_prefix,
                new_prefix=prefix,
            )
        )
        parts.append("  </thead>")

    parts.append("  <tbody>")
    parts.append(
        _render_rows(
            table,
            range(table.n_header_rows, len(table.grid)),
            id_prefix=id_prefix,
            old_prefix=old_prefix,
            new_prefix=prefix,
        )
    )
    parts.append("  </tbody>")

    if include_footnotes:
        foot = _render_footnotes(table, prefix)
        if foot:
            parts.append(foot)

    parts.append("</table>")
    return "\n".join(parts)


def _md_escape(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ").strip()


def table_to_markdown(table: Table, *, include_caption: bool = True) -> str:
    """Render ``table`` as GitHub-flavoured Markdown.

    Markdown has no notion of merged cells, so multi-row headers are flattened into
    composite column keys (``"Proportion ... by Method | Sequence analysis"``). Use
    :func:`table_to_html` when layout fidelity matters.
    """
    keys = table.column_keys()
    if not keys:
        return ""

    lines: list[str] = []
    if include_caption and (table.label or table.caption):
        if table.label:
            lines.append(f"**{table.label.strip()}**")
        if table.caption:
            lines.append(f"*{table.caption}*")
        lines.append("")

    lines.append("| " + " | ".join(_md_escape(k) for k in keys) + " |")
    lines.append("|" + "|".join(" --- " for _ in keys) + "|")
    for row in table.to_grid(filled=True):
        lines.append("| " + " | ".join(_md_escape(v) for v in row) + " |")

    if table.footnotes:
        lines.append("")
        for fn in table.footnotes:
            lines.append(f"- **{fn.label}** {fn.text}")
    return "\n".join(lines)


def footnotes_to_html(footnotes: Iterable[Footnote], *, prefix: str = "fn") -> str:
    """Render a standalone footnote list to ``<ul>`` (for tables rendered without tfoot)."""
    items: list[str] = []
    for fn in footnotes:
        body = fn.html or _attr_text(fn.text)
        if not fn.label:
            items.append(f'<li class="gr-fn-legend">{body}</li>')
            continue
        items.append(
            f'<li id="{_attr_text(footnote_anchor(prefix, fn.label))}">'
            f'<span class="fn-label">{_attr_text(fn.label)}</span> {body}</li>'
        )
    return '<ul class="gr-footnotes">' + "".join(items) + "</ul>" if items else ""


DEFAULT_TABLE_CSS = """\
.gr-table { border-collapse: collapse; margin: 1rem 0; font-size: 0.95rem; }
.gr-table caption { text-align: left; font-weight: 600; padding: 0.35rem 0; }
.gr-table th, .gr-table td { border: 1px solid #b9c0c8; padding: 0.35rem 0.55rem;
  vertical-align: top; text-align: left; }
.gr-table thead th { background: #eef2f6; }
.gr-table tbody th { background: #f7f9fb; font-weight: 600; }
.gr-table .fn-mark a { text-decoration: none; }
.gr-table .fn-orphan { color: #b45309; }  /* marker with no matching footnote */
.gr-footnotes { margin: 0.4rem 0 0; padding-left: 1.2rem; font-size: 0.85rem; }
.gr-footnotes .fn-label { font-weight: 600; }
"""


def html_document(body: str, *, title: str = "GeneReviews", css: str | None = None) -> str:
    """Wrap rendered fragments in a minimal standalone HTML document."""
    style = DEFAULT_TABLE_CSS if css is None else css
    return (
        "<!DOCTYPE html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        f"<title>{_attr_text(title)}</title>\n<style>\n{style}</style>\n"
        f"</head>\n<body>\n{body}\n</body>\n</html>\n"
    )
