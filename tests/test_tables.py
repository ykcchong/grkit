"""Tests for :mod:`grkit.tables` -- span expansion, headers, footnotes, rendering."""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest

from grkit.tables import parse_table_wrap


def _wrap(xml: str):
    return ET.fromstring(xml)


# ──────────────────────────────────────────────────────────────────────────────
# Span expansion: the bug that broke the canonical sequencing-vs-CNV table
# ──────────────────────────────────────────────────────────────────────────────
def test_rowspan_fills_gene_so_yield_stays_in_its_column(table_of):
    """A gene spanning two methods must repeat, not shift the yield into ``Method``.

    Positional parsing returned ``["Gene-targeted deletion/duplication analysis",
    "2.5%"]`` for the second row, which silently attributes a copy-number yield to the
    ``Method`` column. This is the regression guard for that bug.
    """
    table = table_of("tst.T.testing")
    grid = table.to_grid(filled=True)
    assert grid == [
        ["TEST1", "Sequence analysis", "92%"],
        ["TEST1", "Gene-targeted deletion/duplication analysis", "2.5%"],
        ["Unknown", "NA", "5.5%"],
    ]


def test_grid_is_rectangular(table_of):
    table = table_of("tst.T.testing")
    assert table.n_cols == 3
    assert all(len(row) == 3 for row in table.to_grid(filled=True))


def test_to_dicts_keys_follow_headers(table_of):
    table = table_of("tst.T.testing")
    rows = table.to_dicts()
    method_key = "Method"
    yield_key = next(k for k in rows[0] if k.startswith("Proportion of Probands"))
    assert rows[0][method_key] == "Sequence analysis"
    assert rows[1][method_key] == "Gene-targeted deletion/duplication analysis"
    assert rows[1][yield_key] == "2.5%"


def test_stray_header_rowspan_does_not_swallow_first_body_row(table_of):
    """GeneReviews emits ``rowspan="2"`` on single-row ``<thead>`` cells.

    Honouring it literally would push the first data row into columns 2-4 and double the
    width of the table. A header cell may only span *header* rows.
    """
    table = table_of("tst.T.testing")
    assert table.n_header_rows == 1
    assert table.n_cols == 3, "header rowspan must not widen the grid"
    assert table.grid[1][0].text == "TEST1"


def test_body_rowspan_is_honoured(table_of):
    """The gene cell still spans two body rows, and both slots reference one cell."""
    table = table_of("tst.T.testing")
    first = table.cell_at(1, 0)
    assert first.rowspan == 2
    assert table.cell_at(2, 0) is first
    assert first.is_origin(1, 0) and not first.is_origin(2, 0)


# ──────────────────────────────────────────────────────────────────────────────
# Multi-row headers (transposed layout)
# ──────────────────────────────────────────────────────────────────────────────
def test_wide_layout_header_path(table_of):
    table = table_of("tst.T.wide")
    assert table.n_header_rows == 3
    assert table.n_cols == 5
    assert table.column_header_path(2) == [
        "Proportion of Pathogenic Variants Identified by Method",
        "Sequence analysis",
    ]
    assert table.column_header_path(3) == [
        "Proportion of Pathogenic Variants Identified by Method",
        "Gene-targeted deletion/duplication analysis",
        "Common deletions",
    ]


def test_wide_layout_rowspan_does_not_duplicate_labels(table_of):
    """``Gene`` spans all three header rows but must appear once in the path."""
    table = table_of("tst.T.wide")
    assert table.column_header_path(2)[0].startswith("Proportion of Pathogenic")
    assert table.column_header_path(0) == ["Gene"]


def test_header_inside_tbody_is_promoted(table_of):
    """A table whose column headings are ``<th>`` in ``<tbody>`` still gets headers."""
    table = table_of("tst.T.db")
    assert table.n_header_rows == 1
    assert [table.column_leaf_header(c) for c in range(table.n_cols)] == [
        "Gene",
        "Chromosome Locus",
        "Protein",
    ]


# ──────────────────────────────────────────────────────────────────────────────
# Footnotes
# ──────────────────────────────────────────────────────────────────────────────
def test_footnotes_keep_labels_and_text(table_of):
    table = table_of("tst.T.testing")
    labels = [fn.label for fn in table.footnotes]
    assert labels == ["", "1", "2", "3", "4", "5"], "unlabelled legend must not be numbered"
    assert table.footnote_map()["3"].text == "Sequence analysis detects SNVs."


def test_cells_record_their_footnote_markers(table_of):
    table = table_of("tst.T.testing")
    # "2.5% 4, 5" -> markers stripped from text, captured as marks.
    yield_cell = table.cell_at(2, 2)
    assert yield_cell.text == "2.5%"
    assert set(yield_cell.marks) == {"4", "5"}


def test_unit_exponent_is_not_treated_as_a_footnote(table_of):
    """``m<sup>2</sup>`` is attached to a word, so it is a superscript, not a marker."""
    table = table_of("tst.T.units")
    dose = table.cell_at(1, 1)
    assert dose.text == "5 mg/m2"
    assert dose.marks == ()


def test_matrix_exponent_survives_on_a_cell_that_is_only_a_number(table_of):
    """A bare-number cell must not be emptied by marker stripping."""
    table = table_of("tst.T.testing")
    assert table.cell_at(3, 1).text == "NA"


# ──────────────────────────────────────────────────────────────────────────────
# Rendering
# ──────────────────────────────────────────────────────────────────────────────
def test_to_html_is_layout_faithful(table_of):
    html = table_of("tst.T.wide").to_html()
    assert 'rowspan="3"' in html  # Gene spans the header
    assert 'colspan="3"' in html  # "... Identified by Method" group
    assert 'colspan="2"' in html  # del/dup group
    assert 'scope="colgroup"' in html
    assert "<thead>" in html and "<tbody>" in html


def test_to_html_renders_emphasis_once_and_in_order(table_of):
    """Guards against emitting a cell's text before its tag (``TEST1<i></i>``)."""
    html = table_of("tst.T.testing").to_html()
    assert "<i>TEST1</i>" in html
    assert "<i></i>" not in html
    assert "TEST1<i>" not in html


def test_to_html_renders_line_breaks_and_escapes(table_of):
    wide = table_of("tst.T.wide")
    html = wide.to_html()
    assert "Attributed to Pathogenic<br>" in html
    assert "&lt;5%" in html, "a literal < in a cell must be escaped"
    assert wide.cell_at(3, 4).text == "<5%"


def test_to_html_links_markers_to_footnotes(table_of):
    html = table_of("tst.T.testing").to_html()
    assert 'class="fn-ref"' in html
    assert 'href="#fn-tst.T.testing-4"' in html
    assert '<li id="fn-tst.T.testing-4"' in html
    assert "gr-fn-legend" in html, "the unlabelled legend note is still rendered"


def test_to_html_ids_are_unique_and_prefixed(table_of):
    html = table_of("tst.T.testing").to_html(id_prefix="doc1-")
    assert 'id="doc1-tst.T.testing"' in html
    assert "#fn-doc1-tst.T.testing-4" in html


def test_to_markdown_flattens_multi_level_headers(table_of):
    markdown = table_of("tst.T.wide").to_markdown()
    # Pipes inside composite keys are escaped for Markdown.
    assert (
        "Proportion of Pathogenic Variants Identified by Method \\| Sequence analysis"
        in markdown
    )
    assert "| --- " in markdown


def test_table_to_record_shape(table_of):
    from grkit.tables import table_to_record

    record = table_to_record(table_of("tst.T.testing"))
    assert record["n_cols"] == 3
    assert record["n_header_rows"] == 1
    assert record["rows"][1]["Method"] == "Gene-targeted deletion/duplication analysis"
    assert "html" in record
    assert record["footnotes"][1]["label"] == "1"


def test_table_without_table_element_is_safe():
    tw = _wrap('<table-wrap id="x"><caption><p>Empty</p></caption></table-wrap>')
    table = parse_table_wrap(tw)
    assert table.n_cols == 0
    assert table.caption == "Empty"
    assert table.to_html().startswith("<table")


@pytest.mark.parametrize("bad,expected", [("0", 1), ("-2", 1), ("abc", 1), ("", 1), ("3", 3)])
def test_span_attributes_are_sanitised(bad, expected):
    tw = _wrap(
        "<table-wrap><table><thead><tr><th>H</th></tr></thead><tbody>"
        f'<tr><td rowspan="{bad}" colspan="{bad}">v</td></tr>'
        "</tbody></table></table-wrap>"
    )
    table = parse_table_wrap(tw)
    assert table.cell_at(1, 0).rowspan == expected
    assert table.cell_at(1, 0).colspan == expected
