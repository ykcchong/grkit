"""Extract "testing strategy vs yield" from GeneReviews tables, preserving the source.

Almost every GeneReviews chapter carries a table titled *Molecular Genetic Testing Used in
...* whose rows pair a test **method** with the **yield** of that method -- classically
``Sequence analysis`` next to ``Gene-targeted deletion/duplication analysis``, i.e.
sequencing versus copy-number testing. The published layout is not stable, and this module
handles both shapes found in the corpus:

**long** (~696 tables)
    One header row: ``Gene | Method | Proportion ... Detectable by Method``. Each body row
    is one (gene, method, yield) triple. A gene spanning several methods relies on
    ``rowspan``, which :mod:`grkit.tables` already expands correctly.

**wide** (~172 tables)
    A multi-row header in which *methods are columns*::

        Gene | Proportion attributed to gene | Proportion of Variants Identified by Method
                                            | Sequence analysis | Gene-targeted del/dup analysis
                                                                 | Common deletions
                                                                 | Other deletions

    Each body row is one gene; every method column contributes a (gene, method, yield)
    triple, with nested sub-methods (``Common deletions``) kept in their column path.

Yields are messy and are normalized rather than assumed to be percentages: ``100%``,
``<1%``, ``~85%``, ``>98%``, ``5-10%``, ``15/15``, ``NA``, and ``Unknown, none reported``
all occur, and a cell may carry footnote markers (``100% 4, 5``) that are stripped by
:mod:`grkit.inline` before parsing.

The original :class:`~grkit.model.Table` is always retained on the result, so the source
table can still be rendered faithfully as HTML via :meth:`TestingYieldTable.to_html`.
"""

from __future__ import annotations

import csv
import io
import json
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from typing import Any

from .model import Table
from .reader import Chapter

# ──────────────────────────────────────────────────────────────────────────────
# Header vocabulary
# ──────────────────────────────────────────────────────────────────────────────
_METHOD_LEAF_RE = re.compile(r"^(?:test\s+)?method$")
_METHOD_GROUP_RE = re.compile(
    r"identified by (?:test |molecular |genetic )?method"
    r"|detectable by (?:test |molecular |genetic )?method"
    r"|identified by (?:molecular|genetic) testing"
    r"|detection frequency"
    r"|sensitivity$"
)
_YIELD_RE = re.compile(
    r"proportion|detectable|identified|frequency|sensitivity|detection rate|%"
)
_GENE_RE = re.compile(r"^gene\b|^gene$")
_ATTRIBUTION_RE = re.compile(r"attributed to|attributed|contribution|share of")

# Caption variants seen in the corpus: "Molecular Genetic Testing Used in X", "Summary of
# Molecular Genetic Testing Used in X", "Genetic/Genomic Testing Used in X", "Cytogenetic
# and Molecular Genetic Testing Used in X", "X: Molecular Genetic Testing", and the
# "Used for"/"Used to Detect" phrasings.
_CAPTION_RE = re.compile(
    r"(?:molecular\s+genetic|molecular|genetic|genomic|cytogenetic"
    r"(?:\s+and\s+molecular\s+genetic)?)\s+testing"
    r"|molecular\s+genetics?\s*$"
    r"|testing\s+strategy",
    re.I,
)

# Method classification. Order matters: "Gene-targeted deletion/duplication analysis"
# must classify as copy-number, not as "targeted".
_KIND_RULES: tuple[tuple[str, re.Pattern[str]], ...] = (
    (
        "deletion_duplication",
        re.compile(
            r"deletion|duplication|del\s*/\s*dup|mlpa|copy[\s-]number|\bcnv\b|exon[\s-]level"
            r"|intragenic|monosomy|mosaicism",
            re.I,
        ),
    ),
    ("methylation", re.compile(r"methylation", re.I)),
    ("repeat_expansion", re.compile(r"repeat|expansion|trinucleotide|sizing", re.I)),
    # Checked before "sequence" so "Exome sequencing" is not read as single-gene sequencing.
    ("panel_or_genome", re.compile(r"panel|exome|genome|multigene", re.I)),
    ("sequence", re.compile(r"sequence analysis|sequencing|sanger|sequence\s*$", re.I)),
    (
        "cytogenetic",
        re.compile(r"microarray|\bcma\b|karyotype|fish|cytogenetic|chromosom", re.I),
    ),
    ("targeted", re.compile(r"targeted", re.I)),
)

_SEQUENCING_KINDS = frozenset({"sequence"})
_COPY_NUMBER_KINDS = frozenset({"deletion_duplication"})


def _norm(text: str) -> str:
    """Lowercase, collapse whitespace and non-breaking spaces."""
    return re.sub(r"\s+", " ", (text or "").replace("\xa0", " ")).strip().lower()


@dataclass(frozen=True)
class MethodKind:
    """Normalized classification of a testing method."""

    kind: str
    is_sequencing: bool = False
    is_copy_number: bool = False

    @property
    def is_other(self) -> bool:
        return not self.is_sequencing and not self.is_copy_number


def classify_method(label: str) -> MethodKind:
    """Classify a method label into a small vocabulary.

    Returns a :class:`MethodKind` with ``kind`` in ``sequence``, ``deletion_duplication``,
    ``cytogenetic``, ``methylation``, ``repeat_expansion``, ``panel_or_genome``,
    ``targeted``, or ``other``, plus the two flags the sequencing-vs-CNV question needs.
    """
    text = _norm(label)
    if not text:
        return MethodKind("other")
    for kind, pattern in _KIND_RULES:
        if pattern.search(text):
            return MethodKind(
                kind,
                is_sequencing=kind in _SEQUENCING_KINDS,
                is_copy_number=kind in _COPY_NUMBER_KINDS,
            )
    return MethodKind("other")


# ──────────────────────────────────────────────────────────────────────────────
# Yield parsing
# ──────────────────────────────────────────────────────────────────────────────
_QUALIFIER_MAP = {
    "≥": ">=",
    "≤": "<=",
    ">": ">",
    "<": "<",
    "~": "~",
    "≈": "~",
    "about ": "~",
    "up to ": "<=",
    ">=": ">=",
    "<=": "<=",
}
_DASHES = "\u2010\u2011\u2012\u2013\u2014\u2212"

# Words that state a bound or an approximation. They are rewritten to symbols before the
# numeric patterns run, so "Nearly 100%", "Up to 100%" and "Estimated >95%" all parse.
_WORD_QUALIFIERS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"\b(?:greater than|more than|over|in excess of|exceeding)\b", re.I), ">"),
    (re.compile(r"\b(?:at least|no less than|minimum of)\b", re.I), ">="),
    (re.compile(r"\b(?:less than|under|fewer than|below)\b", re.I), "<"),
    (re.compile(r"\b(?:up to|no more than|at most|maximum of)\b", re.I), "<="),
    (
        re.compile(r"\b(?:nearly|almost|approximately|approx\.?|about|around|roughly)\b", re.I),
        "~",
    ),
    (re.compile(r"\b(?:estimated|est\.?)\b", re.I), ""),
)

_QUALIFIER_CLASS = r"[<>~]=?"

# Ranges appear in two shapes in the corpus: "70%-80%" and "5-15%". The second has only a
# trailing percent sign, so it is normalised to the first shape before matching.
_RANGE_NORMALIZE_RE = re.compile(
    rf"(?P<a>\d+(?:\.\d+)?)\s*(?P<dash>[{_DASHES}\-])\s*(?P<b>\d+(?:\.\d+)?)\s*%"
)

# Searched (not anchored) so trailing detail is tolerated: "~97.5% (78/80)".
_PERCENT_SEARCH_RE = re.compile(
    rf"(?P<q>{_QUALIFIER_CLASS})?\s*(?P<v1>\d+(?:\.\d+)?)\s*%"
    rf"(?:\s*[{_DASHES}\-]\s*(?P<q2>{_QUALIFIER_CLASS})?\s*(?P<v2>\d+(?:\.\d+)?)\s*%)?"
)
_FRACTION_SEARCH_RE = re.compile(r"(?P<num>\d+(?:\.\d+)?)\s*/\s*(?P<den>\d+(?:\.\d+)?)")
_OF_RE = re.compile(r"^\s*(?P<num>\d+)\s+of\s+(?P<den>\d+)\b", re.I)
_BARE_NUMBER_RE = re.compile(r"^(?P<q>[<>~]=?)?\s*(?P<v>\d+(?:\.\d+)?)\s*$")
_NA_RE = re.compile(r"^(n/?a|not applicable|none|nil)$", re.I)
_UNKNOWN_RE = re.compile(r"^(unknown|not known|unclear|no data)\b", re.I)
_NONE_REPORTED_RE = re.compile(r"^(none reported|not reported|no .*reported|not observed)\b", re.I)
_ZERO_RE = re.compile(r"^0+\s*%?$|^none detected\b", re.I)
_SEE_FOOTNOTE_RE = re.compile(
    r"^see\s+footnotes?\s+(?P<labels>[0-9][0-9,\s]*?)\s*\.?$", re.I
)

_WORD_NUMBERS = {
    "one": 1.0,
    "two": 2.0,
    "three": 3.0,
    "four": 4.0,
    "five": 5.0,
    "six": 6.0,
    "seven": 7.0,
    "eight": 8.0,
    "nine": 9.0,
    "ten": 10.0,
}
_COUNT_RE = re.compile(
    r"^(?P<n>\d+|one|two|three|four|five|six|seven|eight|nine|ten)\s+"
    r"(?P<unit>reported|individuals?|persons?|people|families|family|pedigrees?|"
    r"cases?|patients?|probands?|variants?|alleles?|report)\b",
    re.I,
)


def _qualifier(raw: str | None) -> str:
    if not raw:
        return ""
    return _QUALIFIER_MAP.get(raw.strip(), raw.strip())


def _apply_word_qualifiers(text: str) -> str:
    for pattern, symbol in _WORD_QUALIFIERS:
        text = pattern.sub(symbol, text)
    text = re.sub(r"\s+", " ", text).strip()
    # "5-15%" -> "5%-15%", so one pattern covers both range spellings.
    return _RANGE_NORMALIZE_RE.sub(
        lambda m: f"{m.group('a')}%{m.group('dash')}{m.group('b')}%", text
    )


@dataclass
class YieldObservation:
    """A parsed yield value.

    ``kind`` is one of:

    ``percent``
        A single percentage, range, or one-sided bound.
    ``fraction``
        A stated ratio (``15/15``) converted to a percentage.
    ``multi_percent``
        Several percentages in one cell (e.g. sporadic vs familial rates); ``raw`` holds
        the full text and the first value is reported.
    ``count``
        A reported number of individuals/families rather than a rate.
    ``see_footnote``
        The cell defers to a table footnote.
    ``not_applicable`` / ``unknown`` / ``none_reported`` / ``empty`` / ``text``
        Non-quantitative outcomes, kept distinct so callers can filter rather than guess.

    ``low``/``high`` bound the value when the source states a limit (``<1%`` -> high=1,
    ``>99%`` -> low=99); ``point`` is set only for an unambiguous single number.
    """

    raw: str = ""
    kind: str = "empty"
    low: float | None = None
    high: float | None = None
    point: float | None = None
    qualifier: str = ""
    numerator: float | None = None
    denominator: float | None = None
    count: float | None = None
    note: str = ""

    @property
    def is_numeric(self) -> bool:
        return self.point is not None or self.low is not None or self.high is not None

    @property
    def percent(self) -> float | None:
        """The unambiguous percent value, or ``None`` for bounds and non-numeric yields."""
        return self.point

    @property
    def midpoint(self) -> float | None:
        """A representative percent for sorting/comparison.

        A range averages; a one-sided bound uses the bound itself (``<1%`` -> 1, ``>98%`` ->
        98), which is the useful approximation when comparing methods. Non-numeric yields
        return ``None``.
        """
        if self.point is not None:
            return self.point
        if self.low is not None and self.high is not None:
            return (self.low + self.high) / 2
        if self.low is not None:
            return self.low
        if self.high is not None:
            return self.high
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "raw": self.raw,
            "kind": self.kind,
            "percent": self.percent,
            "low": self.low,
            "high": self.high,
            "qualifier": self.qualifier,
            "numerator": self.numerator,
            "denominator": self.denominator,
            "count": self.count,
            "note": self.note,
        }


def _fraction_parts(text: str) -> tuple[float | None, float | None]:
    match = _FRACTION_SEARCH_RE.search(text)
    if match:
        return float(match.group("num")), float(match.group("den"))
    match = _OF_RE.match(text)
    if match:
        return float(match.group("num")), float(match.group("den"))
    return None, None


def _percent_from_fraction(num: float | None, den: float | None) -> float | None:
    if num is None or den in (None, 0):
        return None
    return num / den * 100.0


def parse_yield(raw: str, *, assume_percent: bool = True) -> YieldObservation:
    """Parse a yield cell into a :class:`YieldObservation`.

    ``assume_percent`` treats a bare number as a percentage, which is right for a column
    headed "Proportion ... Detectable by Method"; set it to ``False`` to keep bare numbers
    as raw counts.
    """
    text = re.sub(r"\s+", " ", (raw or "").replace("\xa0", " ")).strip()
    if not text:
        return YieldObservation(raw="", kind="empty")

    if _NA_RE.match(text):
        return YieldObservation(raw=text, kind="not_applicable")
    if _ZERO_RE.match(text):
        return YieldObservation(raw=text, kind="percent", low=0.0, high=0.0, point=0.0)
    if _NONE_REPORTED_RE.match(text):
        return YieldObservation(raw=text, kind="none_reported", low=0.0, high=0.0, point=0.0)
    if _UNKNOWN_RE.match(text):
        return YieldObservation(raw=text, kind="unknown")

    footnote = _SEE_FOOTNOTE_RE.match(text)
    if footnote:
        return YieldObservation(
            raw=text, kind="see_footnote", note=footnote.group("labels").strip()
        )

    normalized = _apply_word_qualifiers(text)
    percents = list(_PERCENT_SEARCH_RE.finditer(normalized))

    if percents:
        num, den = _fraction_parts(text)
        first = percents[0]
        qualifier = _qualifier(first.group("q"))
        value = float(first.group("v1"))

        if len(percents) > 1:
            # Several rates in one cell: report the first and keep the text for review.
            return YieldObservation(
                raw=text,
                kind="multi_percent",
                low=value,
                high=value,
                qualifier=qualifier,
                numerator=num,
                denominator=den,
                note=text,
            )

        if first.group("v2") is not None:
            second = float(first.group("v2"))
            return YieldObservation(
                raw=text,
                kind="percent",
                low=min(value, second),
                high=max(value, second),
                qualifier=qualifier,
                numerator=num,
                denominator=den,
            )

        if qualifier in (">", ">="):
            return YieldObservation(
                raw=text, kind="percent", low=value, qualifier=qualifier,
                numerator=num, denominator=den,
            )
        if qualifier in ("<", "<="):
            return YieldObservation(
                raw=text, kind="percent", high=value, qualifier=qualifier,
                numerator=num, denominator=den,
            )
        return YieldObservation(
            raw=text,
            kind="percent",
            low=value,
            high=value,
            point=value,
            qualifier=qualifier,
            numerator=num,
            denominator=den,
        )

    num, den = _fraction_parts(text)
    percent = _percent_from_fraction(num, den)
    if percent is not None:
        return YieldObservation(
            raw=text,
            kind="fraction",
            low=percent,
            high=percent,
            point=percent,
            numerator=num,
            denominator=den,
        )

    count = _COUNT_RE.match(normalized)
    if count:
        token = count.group("n").lower()
        value = _WORD_NUMBERS.get(token, None)
        if value is None:
            value = float(token)
        return YieldObservation(
            raw=text,
            kind="count",
            count=value,
            note=count.group("unit").lower(),
        )

    match = _BARE_NUMBER_RE.match(normalized)
    if match:
        value = float(match.group("v"))
        qualifier = _qualifier(match.group("q"))
        if not assume_percent:
            return YieldObservation(raw=text, kind="count", count=value, point=None)
        if qualifier in (">", ">="):
            return YieldObservation(raw=text, kind="percent", low=value, qualifier=qualifier)
        if qualifier in ("<", "<="):
            return YieldObservation(raw=text, kind="percent", high=value, qualifier=qualifier)
        return YieldObservation(
            raw=text, kind="percent", low=value, high=value, point=value, qualifier=qualifier
        )

    return YieldObservation(raw=text, kind="text")


# ──────────────────────────────────────────────────────────────────────────────
# Result objects
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class YieldRow:
    """One (gene, method, yield) observation."""

    gene: str = ""
    method: str = ""
    method_raw: str = ""
    method_kind: str = "other"
    method_path: list[str] = field(default_factory=list)
    is_sequencing: bool = False
    is_copy_number: bool = False
    yield_raw: str = ""
    observation: YieldObservation = field(default_factory=YieldObservation)
    attribution: str = ""
    detail: str = ""
    marks: tuple[str, ...] = ()
    yield_column: str = ""
    """Leaf header of the column the yield came from (e.g. ``Proband``)."""

    @property
    def yield_percent(self) -> float | None:
        return self.observation.percent

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "method": self.method,
            "method_raw": self.method_raw,
            "method_kind": self.method_kind,
            "is_sequencing": self.is_sequencing,
            "is_copy_number": self.is_copy_number,
            "yield_raw": self.yield_raw,
            "yield_percent": self.yield_percent,
            "yield_midpoint": self.observation.midpoint,
            "yield_low": self.observation.low,
            "yield_high": self.observation.high,
            "yield_qualifier": self.observation.qualifier,
            "yield_status": self.observation.kind,
            "attribution": self.attribution,
            "detail": self.detail,
            "yield_column": self.yield_column,
            "footnote_marks": list(self.marks),
        }


@dataclass
class TestingYieldTable:
    """A testing-strategy-vs-yield table, with the original table preserved."""

    chapter_id: str = ""
    title: str = ""
    table_id: str = ""
    label: str = ""
    caption: str = ""
    section: str = ""
    layout: str = ""
    rows: list[YieldRow] = field(default_factory=list)
    table: Table | None = None
    caption_match: bool = False
    confidence: float = 0.0

    @property
    def genes(self) -> list[str]:
        out: list[str] = []
        for row in self.rows:
            if row.gene and row.gene not in out:
                out.append(row.gene)
        return out

    @property
    def methods(self) -> list[str]:
        out: list[str] = []
        for row in self.rows:
            if row.method and row.method not in out:
                out.append(row.method)
        return out

    @property
    def has_sequencing(self) -> bool:
        return any(r.is_sequencing for r in self.rows)

    @property
    def has_copy_number(self) -> bool:
        return any(r.is_copy_number for r in self.rows)

    def rows_for(self, gene: str) -> list[YieldRow]:
        low = gene.lower()
        return [r for r in self.rows if r.gene.lower() == low]

    def to_html(self, **kwargs: Any) -> str:
        """Render the **original** table faithfully (real rowspan/colspan, footnotes)."""
        if self.table is None:
            return ""
        return self.table.to_html(**kwargs)

    def comparison(self) -> list[dict[str, Any]]:
        """Per-gene sequencing-versus-copy-number summary.

        This is the "testing strategy vs yield" view: for each gene, the sequencing yield
        and the copy-number (deletion/duplication) yield side by side, plus any other
        methods as free-form entries.
        """
        genes: list[str] = []
        for row in self.rows:
            if row.gene not in genes:
                genes.append(row.gene)

        out: list[dict[str, Any]] = []
        for gene in genes:
            rows = self.rows_for(gene)
            record: dict[str, Any] = {"gene": gene}
            sequencing = [r for r in rows if r.is_sequencing]
            copy_number = [r for r in rows if r.is_copy_number]
            other = [r for r in rows if not r.is_sequencing and not r.is_copy_number]
            # Nested sub-methods (e.g. "Common deletions" with an empty cell alongside
            # "Other deletions" with a value) mean the first match may be blank; prefer the
            # first row that actually states a yield.
            seq_row = _first_with_value(sequencing)
            cnv_row = _first_with_value(copy_number)
            record["sequencing"] = seq_row.yield_raw if seq_row else ""
            record["sequencing_percent"] = seq_row.observation.midpoint if seq_row else None
            record["copy_number"] = cnv_row.yield_raw if cnv_row else ""
            record["copy_number_percent"] = cnv_row.observation.midpoint if cnv_row else None
            record["sequencing_methods"] = [
                {"method": r.method, "yield": r.yield_raw, "percent": r.yield_percent}
                for r in sequencing
            ]
            record["copy_number_methods"] = [
                {"method": r.method, "yield": r.yield_raw, "percent": r.yield_percent}
                for r in copy_number
            ]
            record["other_methods"] = [
                {"method": r.method, "yield": r.yield_raw, "column": r.yield_column}
                for r in other
            ]
            record["attribution"] = next((r.attribution for r in rows if r.attribution), "")
            record["methods"] = [
                {
                    "method": r.method,
                    "kind": r.method_kind,
                    "yield": r.yield_raw,
                    "percent": r.yield_percent,
                }
                for r in rows
            ]
            out.append(record)
        return out

    def comparison_to_html(self) -> str:
        """Render the normalized sequencing-vs-copy-number view as an HTML table."""
        from .render import _attr_text

        rows = self.comparison()
        head = (
            "<tr><th scope=\"col\">Gene</th><th scope=\"col\">Sequencing yield</th>"
            "<th scope=\"col\">Copy-number yield</th><th scope=\"col\">Other methods</th></tr>"
        )
        body: list[str] = []
        for record in rows:
            # A method measured in two populations (proband / at-risk relatives) would
            # otherwise render as an indistinguishable "CMA: 100%; CMA: 100%".
            counts: dict[str, int] = {}
            for entry in record["other_methods"]:
                counts[entry["method"]] = counts.get(entry["method"], 0) + 1
            parts: list[str] = []
            for entry in record["other_methods"]:
                label = entry["method"]
                if counts[label] > 1 and entry.get("column"):
                    label = f"{label} ({entry['column']})"
                parts.append(f"{label}: {entry['yield']}")
            other = "; ".join(parts)
            body.append(
                "<tr>"
                f"<th scope=\"row\">{_attr_text(record['gene'])}</th>"
                f"<td>{_attr_text(record['sequencing'])}</td>"
                f"<td>{_attr_text(record['copy_number'])}</td>"
                f"<td>{_attr_text(other)}</td>"
                "</tr>"
            )
        caption = (
            '<caption>Strategy vs yield '
            f'<span class="gr-source">{_attr_text(self.chapter_id)}'
            f' :: {_attr_text(self.table_id)}</span></caption>'
        )
        return (
            '<table class="gr-table gr-strategy">\n'
            f"  {caption}\n  <thead>\n    {head}\n  </thead>\n  <tbody>\n    "
            + "\n    ".join(body)
            + "\n  </tbody>\n</table>"
        )

    def to_dict(self, *, include_source_html: bool = False) -> dict[str, Any]:
        record: dict[str, Any] = {
            "chapter_id": self.chapter_id,
            "title": self.title,
            "table_id": self.table_id,
            "label": self.label,
            "caption": self.caption,
            "section": self.section,
            "layout": self.layout,
            "caption_match": self.caption_match,
            "confidence": self.confidence,
            "genes": self.genes,
            "methods": self.methods,
            "has_sequencing": self.has_sequencing,
            "has_copy_number": self.has_copy_number,
            "rows": [r.to_dict() for r in self.rows],
            "strategy_vs_yield": self.comparison(),
        }
        if include_source_html and self.table is not None:
            record["source_table_html"] = self.to_html()
        return record

    def to_records(self) -> list[dict[str, Any]]:
        """Flat per-observation records, convenient for CSV/DataFrame ingestion."""
        out: list[dict[str, Any]] = []
        for row in self.rows:
            record = {
                "chapter_id": self.chapter_id,
                "title": self.title,
                "table_id": self.table_id,
                "caption": self.caption,
                "layout": self.layout,
                **row.to_dict(),
            }
            out.append(record)
        return out


# ──────────────────────────────────────────────────────────────────────────────
# Extraction
# ──────────────────────────────────────────────────────────────────────────────
def _column_profile(table: Table) -> dict[str, Any]:
    paths = {c: table.column_header_path(c) for c in range(table.n_cols)}
    return {
        "paths": paths,
        "norm_paths": {c: [_norm(p) for p in paths[c]] for c in paths},
        "leaf": {c: _norm(paths[c][-1]) if paths[c] else "" for c in paths},
    }


def _pick_gene_column(profile: dict[str, Any], n_cols: int) -> int | None:
    for c in range(n_cols):
        if _GENE_RE.match(profile["leaf"][c]):
            return c
    # Copy-number-only tables are headed "Deletion 1 | Method | Sensitivity".
    for c in range(n_cols):
        leaf = profile["leaf"][c]
        if leaf and re.match(r"^(deletion|duplication)\b", leaf):
            return c
    return None


def _widen_method_label(path: list[str], group_index: int) -> str:
    """Method label from a column path, dropping the group header at ``group_index``."""
    tail = [p for p in path[group_index + 1 :] if p]
    if not tail:
        tail = [p for p in path if p][-1:]
    return " \u2013 ".join(tail)


def _method_group_index(norm_path: list[str]) -> int | None:
    for i, part in enumerate(norm_path):
        if _METHOD_GROUP_RE.search(part):
            return i
    return None


def _looks_like_method_column(norm_path: list[str], leaf: str) -> bool:
    """True when a header names a *genetic testing method* rather than a data descriptor.

    Deliberately strict, because this is the only path where the method is inferred from a
    header instead of an explicit ``Method`` column or ``... by Method`` group header.
    Gene/database tables head their columns ``Gene``, ``Chromosome Locus``, ``Protein``,
    ``HGMD`` and so on; ``Chromosome Locus`` in particular must not be mistaken for a
    cytogenetic test, or every "Genes and Databases" table becomes a false positive.
    """
    if _METHOD_LEAF_RE.match(leaf) or _ATTRIBUTION_RE.search(leaf):
        return False
    if leaf in _NON_METHOD_HEADERS or any(part in _NON_METHOD_HEADERS for part in norm_path):
        return False
    text = " ".join(norm_path) if norm_path else leaf
    if not _STRONG_METHOD_HINT_RE.search(text):
        return False
    # A test method also has to read like one: "... analysis", "sequencing", "microarray".
    return bool(_TEST_ACTION_RE.search(text))


# Column headings that describe a gene/database record, never a testing method.
_NON_METHOD_HEADERS = frozenset(
    {
        "gene",
        "genes",
        "gene(s)",
        "gene 1",
        "chromosome locus",
        "locus",
        "protein",
        "locus-specific databases",
        "hgmd",
        "clinvar",
        "omim",
        "comment",
        "comments",
        "notes",
        "feature",
        "feature(s)",
        "disorder",
        "moi",
    }
)

_TEST_ACTION_RE = re.compile(
    r"\banalysis\b|\btesting\b|sequencing|\bassay\b|microarray|karyotype|\bfish\b|\bmlpa\b"
    r"|\bpanel\b|screening|cytogenetic",
    re.I,
)


_STRONG_METHOD_HINT_RE = re.compile(
    r"sequence analysis|sequencing|sanger|deletion|duplication|del\s*/\s*dup|mlpa"
    r"|copy[\s-]number|\bcnv\b|microarray|\bcma\b|karyotype|\bfish\b|cytogenetic|methylation"
    r"|repeat expansion|exome|genome|multigene panel",
    re.I,
)


def _cell_text(table: Table, row: int, col: int | None) -> str:
    if col is None:
        return ""
    cell = table.cell_at(row, col)
    return cell.text if cell is not None else ""


def _cell_marks(table: Table, row: int, col: int | None) -> tuple[str, ...]:
    if col is None:
        return ()
    cell = table.cell_at(row, col)
    return cell.marks if cell is not None else ()


def _first_with_value(rows: list[YieldRow]) -> YieldRow | None:
    """Pick the representative row for a method.

    Prefers a quantified yield over a qualitative one, so a nested sub-method that is blank
    or says "Not applicable" does not shadow a sibling that states ``100%``. Falls back to
    the first non-empty row, then the first row.
    """
    for row in rows:
        if row.observation.is_numeric:
            return row
    for row in rows:
        if row.yield_raw.strip():
            return row
    return rows[0] if rows else None


def extract_from_table(table: Table) -> TestingYieldTable | None:
    """Extract a testing-strategy-vs-yield view from one table, or ``None``.

    Returns ``None`` when the table has no recognisable method facet paired with a yield.
    """
    if table.n_cols == 0 or table.n_header_rows == 0 or len(table.body_grid) == 0:
        return None

    profile = _column_profile(table)
    n_cols = table.n_cols
    leaf = profile["leaf"]
    norm_paths = profile["norm_paths"]

    gene_col = _pick_gene_column(profile, n_cols)
    method_cols = [c for c in range(n_cols) if _METHOD_LEAF_RE.match(leaf[c])]
    attribution_cols = [
        c
        for c in range(n_cols)
        if any(_ATTRIBUTION_RE.search(part) for part in norm_paths[c])
        and not any(_METHOD_GROUP_RE.search(part) for part in norm_paths[c])
    ]
    group_cols = [
        c
        for c in range(n_cols)
        if c not in attribution_cols and _method_group_index(norm_paths[c]) is not None
    ]
    plain_yield_cols = [
        c
        for c in range(n_cols)
        if c not in method_cols
        and c != gene_col
        and c not in attribution_cols
        and _YIELD_RE.search(leaf[c])
    ]

    layout: str
    yield_cols: list[int]
    method_labels: dict[int, tuple[str, list[str]]] = {}

    if method_cols:
        layout = "long"
        yield_cols = sorted(set(plain_yield_cols) | set(group_cols))
        yield_cols = [c for c in yield_cols if c not in method_cols and c != gene_col]
        if not yield_cols:
            return None
    elif group_cols:
        layout = "wide"
        yield_cols = [c for c in group_cols if c != gene_col]
        for c in yield_cols:
            idx = _method_group_index(norm_paths[c])
            assert idx is not None
            method_labels[c] = (_widen_method_label(profile["paths"][c], idx), profile["paths"][c])
        if not yield_cols:
            return None
    else:
        # Wide tables with no explicit group header: the columns *are* the methods.
        candidates = [
            c
            for c in range(n_cols)
            if c != gene_col
            and c not in attribution_cols
            and not method_cols
            and _looks_like_method_column(norm_paths[c], leaf[c])
        ]
        if len(candidates) < 1 or gene_col is None:
            return None
        layout = "wide"
        yield_cols = candidates
        for c in yield_cols:
            label = " \u2013 ".join(p for p in profile["paths"][c] if p)
            method_labels[c] = (label, profile["paths"][c])

    if gene_col is None and layout == "wide":
        return None
    if not yield_cols:
        return None

    detail_cols = [
        c
        for c in range(n_cols)
        if c != gene_col
        and c not in yield_cols
        and c not in method_cols
        and c not in attribution_cols
    ]

    body_start = table.n_header_rows
    rows: list[YieldRow] = []
    carry_gene = ""

    for r in range(body_start, len(table.grid)):
        gene = _cell_text(table, r, gene_col).strip() or carry_gene
        if _cell_text(table, r, gene_col).strip():
            carry_gene = gene

        for c in yield_cols:
            raw = _cell_text(table, r, c).strip()
            if layout == "long":
                method_raw = _cell_text(table, r, method_cols[0]).strip()
                method_label = method_raw
                method_path = [method_raw] if method_raw else []
            else:
                method_label, method_path = method_labels.get(c, ("", []))
                method_raw = method_label
            if not raw and not method_label and not gene:
                continue

            kind = classify_method(method_label)
            attribution = ""
            for ac in attribution_cols:
                value = _cell_text(table, r, ac).strip()
                if value:
                    attribution = value
                    break
            detail = " | ".join(
                v for v in (_cell_text(table, r, dc).strip() for dc in detail_cols) if v
            )

            rows.append(
                YieldRow(
                    gene=gene.strip(),
                    method=method_label.strip(),
                    method_raw=method_raw,
                    method_kind=kind.kind,
                    method_path=list(method_path),
                    is_sequencing=kind.is_sequencing,
                    is_copy_number=kind.is_copy_number,
                    yield_raw=raw,
                    observation=parse_yield(raw),
                    attribution=attribution,
                    detail=detail,
                    marks=_cell_marks(table, r, c),
                    yield_column=table.column_leaf_header(c),
                )
            )

    if not rows:
        return None

    numeric = sum(1 for r in rows if r.observation.is_numeric)
    known_method = sum(1 for r in rows if r.method_kind != "other")
    caption_match = bool(_CAPTION_RE.search(table.caption or ""))
    # Confidence answers "is this the right table?", not "are its yields numbers?".
    # A canonical Table 1 whose yields are all "Unknown" is still the right table, so
    # caption and structure carry most of the weight and numeric yield only a little.
    structural = 1.0 if (method_cols or group_cols) else 0.5
    confidence = min(
        1.0,
        0.35 * (1.0 if caption_match else 0.0)
        + 0.30 * structural
        + 0.20 * (known_method / len(rows))
        + 0.15 * (numeric / len(rows)),
    )

    return TestingYieldTable(
        chapter_id=table.chapter_id,
        table_id=table.id,
        label=table.label,
        caption=table.caption,
        section=table.section,
        layout=layout,
        rows=rows,
        table=table,
        caption_match=caption_match,
        confidence=round(confidence, 3),
    )


def is_testing_yield_table(table: Table) -> bool:
    """True when ``table`` pairs a testing method facet with a yield facet."""
    return extract_from_table(table) is not None


def extract_testing_yield(
    chapter: Chapter,
    *,
    require_caption_match: bool = False,
    min_confidence: float = 0.0,
) -> list[TestingYieldTable]:
    """Extract every testing-strategy-vs-yield table from ``chapter``.

    ``require_caption_match`` restricts results to tables whose caption names genetic
    testing, which is the tightest filter for the canonical Table 1. ``min_confidence``
    drops weak structural matches.
    """
    out: list[TestingYieldTable] = []
    for table in chapter.tables:
        result = extract_from_table(table)
        if result is None:
            continue
        result.title = chapter.title
        result.chapter_id = result.chapter_id or chapter.id
        if require_caption_match and not result.caption_match:
            continue
        if result.confidence < min_confidence:
            continue
        out.append(result)
    return out


def extract_all(
    chapters: Iterable[Chapter], **kwargs: Any
) -> Iterator[TestingYieldTable]:
    """Yield testing tables across many chapters."""
    for chapter in chapters:
        yield from extract_testing_yield(chapter, **kwargs)


# ──────────────────────────────────────────────────────────────────────────────
# Export helpers
# ──────────────────────────────────────────────────────────────────────────────
def to_records(tables: Iterable[TestingYieldTable]) -> list[dict[str, Any]]:
    """Flatten many results into per-observation records."""
    out: list[dict[str, Any]] = []
    for table in tables:
        out.extend(table.to_records())
    return out


def write_json(
    tables: Iterable[TestingYieldTable], path: str, *, include_source_html: bool = True
) -> None:
    """Write full results (including preserved table HTML) to JSON."""
    payload = [t.to_dict(include_source_html=include_source_html) for t in tables]
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)


def _csv_value(value: Any) -> Any:
    """Flatten a record value for CSV: lists joined, dicts JSON-encoded, ``None`` blank."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return "; ".join(str(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False)
    return value


def csv_text(records: list[dict[str, Any]]) -> str:
    """Render flat records as CSV text (shared by :func:`write_csv` and the CLI)."""
    if not records:
        return ""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(records[0].keys()))
    writer.writeheader()
    for record in records:
        writer.writerow({key: _csv_value(value) for key, value in record.items()})
    return buffer.getvalue()


def write_csv(tables: Iterable[TestingYieldTable], path: str) -> None:
    """Write flat per-observation records to CSV."""
    text = csv_text(to_records(tables))
    if not text:
        return
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write(text)
