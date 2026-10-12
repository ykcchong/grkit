"""Rich inline content: lossless conversion of NXML inline markup to text and HTML.

GeneReviews table cells are *not* plain strings. A single cell can contain gene names
marked up with ``<italic>``, footnote markers in ``<sup>``, explicit line breaks
(``<break/>``), hyperlinks (``<ext-link>``/``<xref>``/``<related-object>``), and several
block children (``<p>``/``<list-item>``) that must not simply be concatenated.

Naive ``"".join(el.itertext())`` -- what a flat reader does -- loses italics, drops the
distinction between a footnote marker and a unit exponent, and runs sibling blocks
together (``"...syndromeCongenital adrenal hyperplasia..."``). This module fixes that.

Two representations are produced for every element:

``text``
    Whitespace-collapsed plain text with footnote markers removed, safe to parse and
    compare (e.g. a yield of ``"92%"``, not ``"92% 4"``).
``html``
    Escaped inline HTML retaining emphasis, ``<br>``, and links, with footnote markers
    rendered as anchors back to the table's footnote list.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from dataclasses import dataclass, field

# Elements that introduce a visual break between runs of text. Sibling blocks inside one
# cell must be separated, otherwise their text silently concatenates.
_BLOCK_TAGS = frozenset(
    {
        "p",
        "list",
        "list-item",
        "def-list",
        "def-item",
        "term",
        "def",
        "title",
        "label",
        "caption",
        "fn",
    }
)

# Inline elements mapped straight onto their HTML equivalent.
_TAG_TO_HTML = {
    "italic": "i",
    "bold": "b",
    "underline": "u",
    "strike": "s",
    "monospace": "code",
    "sup": "sup",
    "sub": "sub",
    "sc": "span",
    "styled-content": "span",
    "named-content": "span",
}

_LINK_TAGS = frozenset({"ext-link", "uri", "related-object", "related-article"})

# A footnote marker is a short, purely numeric/alphanumeric superscript detached from the
# preceding word. Unit exponents (``m<sup>2</sup>``) are attached, so requiring whitespace
# before the ``<sup>`` separates the two cases.
_MARKER_RE = re.compile(r"^[0-9]{1,2}(?:\s*[,;]\s*[0-9]{1,2})*$|^[a-z]$")
# Collapse every Unicode whitespace, including the &#x000a0; non-breaking spaces that
# GeneReviews uses liberally between a value and its footnote marker.
_WS_RE = re.compile(r"\s+")
_XLINK_HREF = "{http://www.w3.org/1999/xlink}href"


def localname(tag: object) -> str:
    """Return an element's tag without its ``{namespace}`` prefix."""
    if not isinstance(tag, str):
        return ""
    return tag.rsplit("}", 1)[-1]


def _attr(el: ET.Element, *names: str) -> str:
    """Return the first present attribute among ``names`` (namespace-insensitive)."""
    for name in names:
        value = el.get(name)
        if value:
            return value
    for key, value in el.attrib.items():
        if localname(key) in names and value:
            return value
    return ""


def split_markers(text: str) -> tuple[str, tuple[str, ...]]:
    """Split a marker string such as ``"4, 5"`` into ``("4", "5")``."""
    return tuple(p.strip() for p in re.split(r"[,;]", text) if p.strip())


@dataclass
class InlineOptions:
    """Controls how inline content is rendered.

    ``footnote_labels`` is the set of footnote labels actually defined for the enclosing
    table. A ``<sup>`` is only treated as a footnote reference when its value is one of
    these labels *and* it is detached from the preceding word; everything else stays a
    plain superscript. This is data-driven, so tables without footnotes never lose text.
    """

    footnote_labels: frozenset[str] = frozenset()
    anchor_prefix: str = "fn"
    resolve_link: object | None = None  # callable(el) -> str href
    extra_classes: dict[str, str] = field(default_factory=dict)


# ──────────────────────────────────────────────────────────────────────────────
# Plain text
# ──────────────────────────────────────────────────────────────────────────────
def element_to_text(el: ET.Element | None) -> str:
    """Whitespace-collapsed plain text of ``el``, with block boundaries respected."""
    if el is None:
        return ""
    out: list[str] = []
    _emit(el, out, "", frozenset(), [], False)
    return _WS_RE.sub(" ", "".join(out)).strip()


# ──────────────────────────────────────────────────────────────────────────────
# HTML
# ──────────────────────────────────────────────────────────────────────────────
def _href_for(el: ET.Element) -> str:
    href = _attr(el, "href", "xlink:href", "uri")
    if href:
        return href
    if localname(el.tag) == "xref":
        rid = _attr(el, "rid")
        if rid:
            return f"#{rid}"
    return ""


def _footnote_anchor(prefix: str, label: str) -> str:
    return f"{prefix}-{label}" if prefix else label


def _tags_for(el: ET.Element, opts: InlineOptions) -> tuple[str, str, bool]:
    """Return ``(open_tag, close_tag, self_contained)`` for ``el``."""
    name = localname(el.tag)
    if name in ("break", "br"):
        return "<br>", "", True
    if name in _TAG_TO_HTML:
        html_name = _TAG_TO_HTML[name]
        attrs = ""
        if name == "sc":
            attrs = ' class="small-caps"'
        elif name == "named-content":
            attrs = ' class="named-content"'
        elif name in opts.extra_classes:
            attrs = f' class="{opts.extra_classes[name]}"'
        return f"<{html_name}{attrs}>", f"</{html_name}>", False
    if name in _LINK_TAGS or name == "xref":
        href = _href_for(el)
        ref_type = _attr(el, "ref-type")
        if href:
            cls = "xref" if name == "xref" else "ext-link"
            if ref_type:
                cls += f" xref-{ref_type}"
            return f'<a class="{cls}" href="{_escape(href, quote=True)}">', "</a>", False
        return '<span class="link">', "</span>", False
    if name in ("p", "term", "def", "list-item", "def-item"):
        tag = "p" if name in ("p", "term", "def") else "div"
        return f"<{tag}>", f"</{tag}>", False
    if name == "list":
        tag = "ol" if _attr(el, "list-type") == "ordered" else "ul"
        return f"<{tag}>", f"</{tag}>", False
    if name == "graphic":
        href = _href_for(el)
        if href:
            return f'<img src="{_escape(href, quote=True)}" alt="">', "", True
        return "", "", True
    return "", "", False


def _render_html(
    el: ET.Element,
    out: list[str],
    opts: InlineOptions,
    *,
    preceding: str,
) -> None:
    """Render ``el`` and its children into ``out``.

    ``preceding`` is the raw text emitted *before* this element. It is essential that it
    excludes the element's own text: a footnote marker is recognised precisely because the
    text before the ``<sup>`` ends in whitespace, which is what separates ``Gene&nbsp;<sup>1
    </sup>`` (a marker) from ``m<sup>2</sup>`` (a unit exponent).
    """
    name = localname(el.tag)

    if opts.footnote_labels and name == "sup":
        plain = element_to_text(el)
        if plain and _is_detached(preceding) and _MARKER_RE.match(plain):
            labels = split_markers(plain)
            if labels and all(lbl in opts.footnote_labels for lbl in labels):
                links = ", ".join(
                    f'<a class="fn-ref" href="#'
                    f'{_footnote_anchor(opts.anchor_prefix, lbl)}">{lbl}</a>'
                    for lbl in labels
                )
                out.append(f'<sup class="fn-mark">{links}</sup>')
                if el.tail:
                    out.append(_escape(el.tail))
                return
            if labels:
                # Orphan marker: the table has footnotes but not this label. Keep it
                # visible (the rendering stays faithful to the published table) while
                # flagging it, since there is no footnote to link to.
                out.append(f'<sup class="fn-orphan">{_escape(plain)}</sup>')
                if el.tail:
                    out.append(_escape(el.tail))
                return

    open_tag, close_tag, self_contained = _tags_for(el, opts)
    out.append(open_tag)
    if not self_contained:
        if el.text:
            out.append(_escape(el.text))
        inner = preceding + (el.text or "")
        for child in el:
            _render_html(child, out, opts, preceding=inner)
            inner += element_to_text(child) + (child.tail or "")
        out.append(close_tag)
    if el.tail:
        out.append(_escape(el.tail))


def _is_detached(preceding: str) -> bool:
    """True when nothing (or only whitespace/punctuation) precedes the element.

    ``m<sup>2</sup>``  -> attached   -> exponent
    ``Gene&nbsp;<sup>1</sup>`` -> detached -> footnote marker
    """
    if not preceding:
        return True
    return preceding[-1].isspace() or preceding[-1] in "%):;,."


def _escape(text: str, *, quote: bool = False) -> str:
    if not text:
        return ""
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    if quote:
        text = text.replace('"', "&quot;")
    return text


def element_to_html(el: ET.Element | None, opts: InlineOptions | None = None) -> str:
    """Return escaped inline HTML for ``el`` (emphasis, breaks, links, footnote anchors)."""
    if el is None:
        return ""
    opts = opts or InlineOptions()
    out: list[str] = []
    if el.text:
        out.append(_escape(el.text))
    preceding = el.text or ""
    for child in el:
        _render_html(child, out, opts, preceding=preceding)
        preceding += element_to_text(child) + (child.tail or "")
    return _WS_RE.sub(" ", "".join(out)).strip()


# ──────────────────────────────────────────────────────────────────────────────
# Marker-aware plain text
# ──────────────────────────────────────────────────────────────────────────────
def _emit(
    node: ET.Element,
    out: list[str],
    preceding: str,
    known: frozenset[str],
    marks: list[str],
    allow_marks: bool,
    orphans: list[str] | None = None,
) -> str:
    """Append ``node``'s text to ``out``, dropping footnote markers.

    Returns the raw text emitted so far, so the caller can keep deciding whether a
    following ``<sup>`` is detached.

    Three cases are distinguished for a detached, marker-shaped ``<sup>``:

    * label is defined in the table -> removed from the text, recorded in ``marks``, and
      linked in the HTML rendering;
    * label is *not* defined (an orphan: GeneReviews sometimes drops a footnote while
      leaving its marker behind, e.g. ``brugada`` defines 1-5, 7, 8 but a cell cites 6) ->
      still removed from the text so parsed values stay clean, but recorded in ``orphans``
      so the gap is visible rather than silently swallowed;
    * not marker-shaped or not detached -> kept verbatim (a unit exponent such as
      ``m<sup>2</sup>``, or a bare number that is the cell's own value).

    ``raw_text`` on the cell always preserves the original, so nothing is lost.
    """
    if node.text:
        out.append(node.text)
        preceding += node.text

    for child in node:
        name = localname(child.tag)
        if allow_marks and name == "sup":
            value = element_to_text(child)
            if value and _is_detached(preceding) and _MARKER_RE.match(value):
                labels = split_markers(value)
                if labels:
                    confirmed = [lbl for lbl in labels if lbl in known]
                    missing = [lbl for lbl in labels if lbl not in known]
                    for lbl in confirmed:
                        if lbl not in marks:
                            marks.append(lbl)
                    if orphans is not None:
                        for lbl in missing:
                            if lbl not in orphans:
                                orphans.append(lbl)
                    if child.tail:
                        out.append(child.tail)
                        preceding += child.tail
                    continue
        if name in ("break", "br"):
            out.append(" ")
            preceding += " "
        elif name in _BLOCK_TAGS:
            out.append(" ")
            preceding += " "
            preceding = _emit(child, out, preceding, known, marks, allow_marks, orphans)
            out.append(" ")
            preceding += " "
        else:
            preceding = _emit(child, out, preceding, known, marks, allow_marks, orphans)
        if child.tail:
            out.append(child.tail)
            preceding += child.tail
    return preceding


def element_to_text_and_marks(
    el: ET.Element | None,
    labels: Iterable[str] = (),
    orphans: list[str] | None = None,
) -> tuple[str, tuple[str, ...]]:
    """Return ``(text, marks)`` for ``el``.

    ``text`` is whitespace-collapsed with footnote markers removed, so a yield reads
    ``"92%"`` rather than ``"92% 4"``. Tags whose labels are not defined in the table are
    removed too, and collected into ``orphans`` when a list is supplied. A cell whose whole
    content is a bare number (``"1"``) is preserved, so score columns are never emptied.
    """
    if el is None:
        return "", ()
    known = frozenset(labels)
    out: list[str] = []
    marks: list[str] = []
    _emit(el, out, "", known, marks, bool(known), orphans)
    return _WS_RE.sub(" ", "".join(out)).strip(), tuple(marks)


def footnote_marks(el: ET.Element | None, labels: Iterable[str]) -> tuple[str, ...]:
    """Return the footnote labels referenced by detached ``<sup>`` markers inside ``el``."""
    return element_to_text_and_marks(el, labels)[1]


def strip_marks(text: str, labels: Iterable[str]) -> str:
    """Remove trailing/embedded footnote markers from already-flattened ``text``.

    Safety net for text that was flattened before marker handling, e.g. ``"100% 4, 5"``.
    Only digits that are known footnote labels are removed, so real values survive.
    """
    known = {lbl for lbl in labels}
    if not known or not text:
        return text
    tokens = text.split()
    kept = [
        t
        for t in tokens
        if not (re.fullmatch(r"[0-9]{1,2}[,;]?", t) and t.rstrip(",;") in known)
    ]
    return " ".join(kept).strip()
