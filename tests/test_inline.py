"""Tests for :mod:`grkit.inline` -- text extraction and inline HTML fidelity."""

from __future__ import annotations

import xml.etree.ElementTree as ET

from grkit.inline import (
    InlineOptions,
    element_to_html,
    element_to_text,
    element_to_text_and_marks,
    strip_marks,
)


def el(xml: str) -> ET.Element:
    return ET.fromstring(xml)


# ──────────────────────────────────────────────────────────────────────────────
# Plain text
# ──────────────────────────────────────────────────────────────────────────────
def test_sibling_blocks_are_separated():
    """Sibling <p> in one cell must not concatenate into one word soup."""
    text = element_to_text(
        el("<caption><p>Antley-Bixler syndrome</p><p>Congenital adrenal hyperplasia</p></caption>")
    )
    assert text == "Antley-Bixler syndrome Congenital adrenal hyperplasia"


def test_break_becomes_a_space_in_text():
    assert element_to_text(el("<p>Alpha<break/>Beta</p>")) == "Alpha Beta"


def test_nested_inline_text_is_joined():
    assert element_to_text(el("<p>the <italic>TEST1</italic> gene</p>")) == "the TEST1 gene"


def test_entities_and_nbsp_are_normalised():
    assert element_to_text(el("<p>a&#x000a0;b &amp; c</p>")) == "a b & c"


# ──────────────────────────────────────────────────────────────────────────────
# Footnote markers
# ──────────────────────────────────────────────────────────────────────────────
def test_detached_known_label_is_removed_and_recorded():
    text, marks = element_to_text_and_marks(el("<td>92%&#x000a0;<sup>4</sup></td>"), {"4"})
    assert text == "92%"
    assert marks == ("4",)


def test_multiple_labels_are_all_recorded():
    text, marks = element_to_text_and_marks(
        el("<td>2.5%&#x000a0;<sup>4,&#x000a0;5</sup></td>"), {"4", "5"}
    )
    assert text == "2.5%"
    assert set(marks) == {"4", "5"}


def test_orphan_marker_is_excluded_from_text_but_recorded():
    """A marker whose footnote is missing still must not pollute the parsed value.

    Real case: the ``brugada`` testing table defines labels 1-5, 7, 8 yet a cell cites 6,
    leaving "~1% 6" in the flattened text. The value is kept clean and the gap recorded.
    """
    orphans: list[str] = []
    text, marks = element_to_text_and_marks(
        el("<td>92%&#x000a0;<sup>9</sup></td>"), {"4"}, orphans
    )
    assert text == "92%"
    assert marks == ()
    assert orphans == ["9"]


def test_unit_exponent_is_not_a_marker():
    """``m<sup>2</sup>`` is attached to a word, so it stays a superscript."""
    text, marks = element_to_text_and_marks(el("<td>5 mg/m<sup>2</sup></td>"), {"2"})
    assert text == "5 mg/m2"
    assert marks == ()


def test_bare_number_cell_is_preserved():
    """A score column containing '1' must not be emptied by marker stripping."""
    text, marks = element_to_text_and_marks(el("<td>1</td>"), {"1", "2"})
    assert text == "1"
    assert marks == ()


def test_marker_at_start_of_cell_is_detected():
    text, marks = element_to_text_and_marks(el("<td><sup>1</sup></td>"), {"1"})
    assert text == ""
    assert marks == ("1",)


def test_orphan_marker_renders_unlinked_and_flagged():
    opts = InlineOptions(footnote_labels=frozenset({"4"}), anchor_prefix="fn-tbl")
    html = element_to_html(el("<td>~1%&#x000a0;<sup>6</sup></td>"), opts)
    assert html == '~1% <sup class="fn-orphan">6</sup>'
    assert "fn-ref" not in html


def test_no_labels_means_no_marker_handling():
    text, marks = element_to_text_and_marks(el("<td>92%&#x000a0;<sup>4</sup></td>"), ())
    assert text == "92% 4"
    assert marks == ()


def test_strip_marks_safety_net():
    assert strip_marks("100% 4, 5", {"4", "5"}) == "100%"
    assert strip_marks("100% 4", {"7"}) == "100% 4"


# ──────────────────────────────────────────────────────────────────────────────
# Inline HTML
# ──────────────────────────────────────────────────────────────────────────────
def test_emphasis_tag_wraps_its_text():
    """Guards the ordering bug that produced ``TEST1<i></i>``.

    The outermost element is rendered inline-only (its own tag is supplied by the table
    renderer), so only the children's tags appear here.
    """
    assert element_to_html(el("<p><italic toggle='yes'>TEST1</italic> gene</p>")) == (
        "<i>TEST1</i> gene"
    )


def test_tail_text_is_not_duplicated():
    assert element_to_html(el("<p>a<italic>b</italic>c<bold>d</bold>e</p>")) == (
        "a<i>b</i>c<b>d</b>e"
    )


def test_break_becomes_br():
    assert element_to_html(el("<td>Alpha<break/>Beta</td>")) == "Alpha<br>Beta"


def test_marker_becomes_an_anchor():
    opts = InlineOptions(footnote_labels=frozenset({"4"}), anchor_prefix="fn-tbl")
    html = element_to_html(el("<td>92%&#x000a0;<sup>4</sup></td>"), opts)
    assert html == '92% <sup class="fn-mark"><a class="fn-ref" href="#fn-tbl-4">4</a></sup>'


def test_marker_with_two_labels_links_both():
    opts = InlineOptions(footnote_labels=frozenset({"4", "5"}), anchor_prefix="fn-tbl")
    html = element_to_html(el("<td>x&#x000a0;<sup>4, 5</sup></td>"), opts)
    assert html.count('class="fn-ref"') == 2
    assert 'href="#fn-tbl-4"' in html and 'href="#fn-tbl-5"' in html


def test_links_are_rendered_and_escaped():
    source = el(
        "<p><ext-link href='https://x.test/a?b=1&amp;c=2'>A &amp; B</ext-link></p>"
    )
    html = element_to_html(source)
    assert 'href="https://x.test/a?b=1&amp;c=2"' in html
    assert "A &amp; B" in html


def test_markup_characters_in_text_are_escaped():
    html = element_to_html(el("<td>&lt;5% &amp; &gt;1%</td>"))
    assert html == "&lt;5% &amp; &gt;1%"


def test_unhandled_inline_elements_are_transparent():
    assert element_to_html(el("<td>a<named-content>x</named-content>b</td>")) == (
        'a<span class="named-content">x</span>b'
    )


def test_none_is_safe():
    assert element_to_text(None) == ""
    assert element_to_html(None) == ""
    assert element_to_text_and_marks(None) == ("", ())
