"""Gene -> testing-coverage tables.

Answers "for gene X, what testing is recommended and what is its expected yield?" by
projecting the corpus onto a single gene:

    >>> from grkit import Corpus, gene_testing_coverage
    >>> with Corpus("gene_NBK1116.tar.gz") as corpus:
    ...     cov = gene_testing_coverage("SCN1A", corpus.iter_chapters())
    >>> cov.summary_rows()
    [{'method': 'Sequence analysis', 'kind': 'sequence', 'entries': 2, ...}]

Three things make this more than a filter over :mod:`grkit.testing`:

**Cell-level gene matching.** A GeneReviews gene cell is not one symbol. It can be
``HBA1 & HBA2``, or a fifteen-symbol mitochondrial run
(``MT-TC MT-TF MT-TH ...``), or ``Unknown``. :func:`cell_genes` splits and normalizes
those, so a query for ``HBA2`` finds the ``HBA1 & HBA2`` row instead of missing it.

**Implicit attribution.** Chapters such as ``gefs`` (SCN1A Seizure Disorders) have a
testing table with no usable Gene column because the chapter is about exactly one gene.
When a chapter indexes a single gene, its rows are attributed to that gene and flagged
``implicit`` so the inference is visible.

**Coverage gaps are reported, not hidden.** An entry can index the gene yet carry no
testing table at all. For SCN5A that is the normal case: three entries index it, but only
``brugada`` has a testing table (the overviews ``dcm-ov`` and ``rws`` do not). Those
entries appear with ``has_table=False`` rather than silently vanishing, because "this
entry discusses the gene but publishes no method/yield data" is itself the answer.
"""

from __future__ import annotations

import re
import statistics
from collections.abc import Iterable
from dataclasses import dataclass, field
from html import escape
from typing import Any

from .reader import Chapter
from .testing import TestingYieldTable, YieldObservation, extract_testing_yield

# Gene cells join symbols with "&", "/", ",", ";" or whitespace, and can be padded or
# wrapped in parentheses. These are never gene symbols.
_NON_GENES = frozenset(
    {
        "",
        "-",
        "na",
        "n/a",
        "none",
        "none reported",
        "not applicable",
        "unknown",
        "various",
        "multiple",
        "several",
        "other",
    }
)
# Individual words that are never a gene symbol but survive the split, e.g. the tail of
# "None reported" or "No variants reported" in a Gene cell.
_GENE_STOPWORDS = frozenset(
    {
        "reported",
        "not",
        "no",
        "and",
        "or",
        "variants",
        "variant",
        "mutations",
        "mutation",
        "detected",
        "see",
        "footnote",
        "unknown",
        "other",
        "others",
    }
)

# Some Gene columns hold a sentence rather than symbols -- the recurrent-deletion entries
# label their rows "~0.8-Mb heterozygous deletion at 1q21.1". Splitting that yields
# "HETEROZYGOUS", "SEQ[GRCH38" and similar, so a token must look like a gene symbol:
# a letter followed by alphanumerics, optionally hyphen/dot separated (MT-TC, HLA-DRB1,
# C11orf65, ATP6V0A2). Chromosome arms ("1Q21.1") and markup fragments ("ISCN:", "SEQ[")
# are rejected by shape.
_GENE_TOKEN_RE = re.compile(r"^[A-Z][A-Z0-9]*(?:[-.][A-Z0-9]+)*$")

# The C11orf65 / C2orf71 family: symbols whose lowercase "orf" is legitimate.
_ORF_LIKE_RE = re.compile(r"^[A-Z]?\d*orf\d+$", re.I)

# Prefixes used by curation bodies and assemblies, never by gene symbols.
_NON_GENE_PREFIXES = ("ISCA", "CLINGEN", "GRCH", "HG19", "HG38")

# Uppercase words that match the shape above but are prose or assay abbreviations.
_PROSE_TOKENS = frozenset(
    {
        "AT",
        "BP",
        "CMA",
        "CNV",
        "CNVS",
        "DEL",
        "DELS",
        "DELETION",
        "DELETIONS",
        "DUP",
        "DUPS",
        "DUPLICATION",
        "DUPLICATIONS",
        "FISH",
        "GENE",
        "GENES",
        "GRCH37",
        "GRCH38",
        "HG19",
        "HG38",
        "HETEROZYGOUS",
        "HOMOZYGOUS",
        "INDEL",
        "ID",
        "IDS",
        "ISCN",
        "KB",
        "MB",
        "MLPA",
        "NGS",
        "PCR",
        "QPCR",
        "SEQ",
        "SNV",
        "SNVS",
        "WES",
        "WGS",
        # Sex-chromosome complements from karyotype notation ("46,XY"), not genes.
        "XX",
        "XXX",
        "XXY",
        "XY",
        "XYY",
    }
)
_GENE_SPLIT_RE = re.compile(r"[,;/&]|\band\b|\s+", re.I)
_SEQ_LABELS = {"sequence": "Sequence analysis"}
_CNV_LABELS = {"deletion_duplication": "Deletion/duplication analysis"}


def _looks_like_symbol(raw: str) -> bool:
    """True when a raw token can be a gene symbol rather than a word of prose.

    Casing is the decisive signal and must be checked *before* uppercasing: GeneReviews
    writes symbols in caps (``HBA1``, ``MT-TC``) while sentences are lowercase
    (``recurrent ... to ... hybrid allele``). The one common exception is the
    ``C11orf65``/``C2orf71`` family, which carries a lowercase ``orf``.
    """
    token = raw.strip(" .,()[]")
    if not token:
        return False
    if _ORF_LIKE_RE.match(token):
        return True
    return not any(ch.islower() for ch in token)


def cell_genes(label: str) -> list[str]:
    """Split a gene cell into normalized, de-duplicated symbols.

    ``"HBA1 & HBA2"`` -> ``["HBA1", "HBA2"]`` (nbsp-safe). Non-symbols are dropped:
    placeholders (``Unknown``, ``NA``, ``-``), prose from sentence-like cells
    (``recurrent ... hybrid allele``), chromosome arms (``1Q21.1``), assay abbreviations
    (``MLPA``) and region ids (``ISCA-37392``).
    """
    if not label:
        return []
    text = re.sub(r"\s+", " ", label.replace("\xa0", " ")).strip()
    if text.strip(" .,;:").lower() in _NON_GENES:
        return []
    out: list[str] = []
    for part in _GENE_SPLIT_RE.split(text):
        if not _looks_like_symbol(part):
            continue
        token = part.strip(" .,()[]").upper()
        if not token or token.lower() in _NON_GENES or token.lower() in _GENE_STOPWORDS:
            continue
        if token in _PROSE_TOKENS or not _GENE_TOKEN_RE.match(token):
            continue
        if token.startswith(_NON_GENE_PREFIXES):
            continue
        if token not in out:
            out.append(token)
    return out


@dataclass
class CoverageRow:
    """One (entry, method, yield) observation for the queried gene."""

    gene: str
    gene_label: str
    chapter_id: str
    chapter_title: str
    table_id: str
    caption: str = ""
    method: str = ""
    method_raw: str = ""
    method_kind: str = "other"
    is_sequencing: bool = False
    is_copy_number: bool = False
    yield_raw: str = ""
    observation: YieldObservation = field(default_factory=YieldObservation)
    attribution: str = ""
    detail: str = ""
    order: int = 0
    """Position within the source table; GeneReviews lists the preferred strategy first."""

    implicit: bool = False
    """True when the gene came from the chapter, not from a Gene column."""

    @property
    def yield_percent(self) -> float | None:
        return self.observation.percent

    @property
    def yield_midpoint(self) -> float | None:
        return self.observation.midpoint

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "gene_label": self.gene_label,
            "chapter_id": self.chapter_id,
            "chapter_title": self.chapter_title,
            "table_id": self.table_id,
            "method": self.method,
            "method_kind": self.method_kind,
            "is_sequencing": self.is_sequencing,
            "is_copy_number": self.is_copy_number,
            "yield_raw": self.yield_raw,
            "yield_percent": self.yield_percent,
            "yield_midpoint": self.yield_midpoint,
            "yield_status": self.observation.kind,
            "attribution": self.attribution,
            "detail": self.detail,
            "order": self.order,
            "implicit": self.implicit,
        }


@dataclass
class CoverageEntry:
    """An entry that either indexes the gene, publishes method/yield data for it, or both."""

    chapter_id: str
    title: str
    indexed: bool = False
    has_table: bool = False
    table_id: str = ""
    layout: str = ""
    confidence: float = 0.0
    row_count: int = 0
    gene_id: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "chapter_id": self.chapter_id,
            "title": self.title,
            "indexed": self.indexed,
            "has_table": self.has_table,
            "table_id": self.table_id,
            "layout": self.layout,
            "confidence": self.confidence,
            "row_count": self.row_count,
            "gene_id": self.gene_id,
        }


def _median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


@dataclass
class GeneTestingCoverage:
    """A gene's testing coverage across the corpus."""

    gene: str = ""
    gene_ids: list[str] = field(default_factory=list)
    aliases: list[str] = field(default_factory=list)
    entries: list[CoverageEntry] = field(default_factory=list)
    rows: list[CoverageRow] = field(default_factory=list)

    # ── shape ────────────────────────────────────────────────────────────────
    @property
    def has_coverage(self) -> bool:
        return bool(self.rows)

    @property
    def entries_with_table(self) -> list[CoverageEntry]:
        return [e for e in self.entries if e.has_table]

    @property
    def entries_indexed_only(self) -> list[CoverageEntry]:
        """Entries that discuss the gene but publish no method/yield table for it."""
        return [e for e in self.entries if e.indexed and not e.has_table]

    @property
    def n_entries(self) -> int:
        return len({r.chapter_id for r in self.rows})

    def sequencing(self) -> list[CoverageRow]:
        return [r for r in self.rows if r.is_sequencing]

    def copy_number(self) -> list[CoverageRow]:
        return [r for r in self.rows if r.is_copy_number]

    def rows_for(self, chapter_id: str) -> list[CoverageRow]:
        return [r for r in self.rows if r.chapter_id == chapter_id]

    def first_line(self) -> CoverageRow | None:
        """The first-listed method row, which is the preferred strategy in GeneReviews."""
        return min(self.rows, key=lambda r: r.order) if self.rows else None

    def method_kinds(self) -> list[str]:
        out: list[str] = []
        for row in self.rows:
            if row.method_kind not in out:
                out.append(row.method_kind)
        return out

    # ── aggregation ──────────────────────────────────────────────────────────
    def summary_rows(self) -> list[dict[str, Any]]:
        """One row per method, with the yield range observed across entries."""
        buckets: dict[str, list[CoverageRow]] = {}
        for row in self.rows:
            buckets.setdefault(row.method, []).append(row)

        out: list[dict[str, Any]] = []
        for method, rows in buckets.items():
            # "None reported" is an explicit zero and joins the numeric range; anything
            # else that is not a rate ("Unknown", "1 family", "See footnote 7.") is kept
            # out of the range and surfaced separately so the range is never read as the
            # whole story.
            numeric = [
                r.yield_midpoint
                for r in rows
                if r.observation.kind in ("percent", "fraction") and r.yield_midpoint is not None
            ]
            qualitative = [
                r.yield_raw
                for r in rows
                if r.observation.kind not in ("percent", "fraction", "none_reported")
            ]
            zero_reported = [
                r.yield_raw for r in rows if r.observation.kind == "none_reported"
            ]
            percents = numeric + [0.0] * len(zero_reported)
            out.append(
                {
                    "method": method,
                    "method_kind": rows[0].method_kind,
                    "is_sequencing": rows[0].is_sequencing,
                    "is_copy_number": rows[0].is_copy_number,
                    "entries": len({r.chapter_id for r in rows}),
                    "observations": len(rows),
                    "yields": [r.yield_raw for r in rows],
                    "min_percent": min(percents) if percents else None,
                    "max_percent": max(percents) if percents else None,
                    "median_percent": _median(percents),
                    "qualitative": qualitative,
                    "first_order": min(r.order for r in rows),
                }
            )
        out.sort(key=lambda r: r["first_order"])
        return out

    def coverage_matrix(self) -> dict[str, dict[str, str]]:
        """``{chapter_id: {method: yield}}`` -- the raw gene x method grid."""
        matrix: dict[str, dict[str, str]] = {}
        for row in self.rows:
            matrix.setdefault(row.chapter_id, {})[row.method] = row.yield_raw
        return matrix

    def sensitivity(self) -> dict[str, Any]:
        """Sequencing versus copy-number sensitivity, the headline coverage question.

        Returns the observed yields and the best numeric estimate per arm, plus which arm
        the source tables list first (GeneReviews puts the preferred strategy first).
        """
        seq, cnv = self.sequencing(), self.copy_number()

        def arm(rows: list[CoverageRow]) -> dict[str, Any]:
            numeric = [
                r.yield_midpoint
                for r in rows
                if r.observation.kind in ("percent", "fraction") and r.yield_midpoint is not None
            ]
            qualitative = [
                r.yield_raw
                for r in rows
                if r.observation.kind not in ("percent", "fraction", "none_reported")
            ]
            return {
                "yields": [r.yield_raw for r in rows],
                "entries": sorted({r.chapter_id for r in rows}),
                "best_percent": max(numeric) if numeric else None,
                "median_percent": _median(numeric),
                "qualitative": qualitative,
                "first_order": min((r.order for r in rows), default=None),
            }

        seq_arm, cnv_arm = arm(seq), arm(cnv)
        preferred = None
        if seq_arm["first_order"] is not None or cnv_arm["first_order"] is not None:
            seq_first = seq_arm["first_order"] if seq_arm["first_order"] is not None else 10**9
            cnv_first = cnv_arm["first_order"] if cnv_arm["first_order"] is not None else 10**9
            if seq_first != cnv_first:
                preferred = "sequence" if seq_first < cnv_first else "deletion_duplication"
        return {"sequencing": seq_arm, "copy_number": cnv_arm, "preferred": preferred}

    # ── serialization ────────────────────────────────────────────────────────
    def to_records(self) -> list[dict[str, Any]]:
        return [row.to_dict() for row in self.rows]

    def to_dict(self, *, include_source_html: bool = False) -> dict[str, Any]:
        record: dict[str, Any] = {
            "gene": self.gene,
            "gene_ids": self.gene_ids,
            "aliases": self.aliases,
            "entries_total": len(self.entries),
            "entries_with_table": len(self.entries_with_table),
            "entries_indexed_only": [e.chapter_id for e in self.entries_indexed_only],
            "has_coverage": self.has_coverage,
            "methods": self.method_kinds(),
            "summary": self.summary_rows(),
            "sensitivity": self.sensitivity(),
            "matrix": self.coverage_matrix(),
            "entries": [e.to_dict() for e in self.entries],
            "rows": self.to_records(),
        }
        if include_source_html:
            record["source_tables_html"] = self.source_html()
        return record

    # ── rendering ────────────────────────────────────────────────────────────
    def source_html(self, results: list[TestingYieldTable] | None = None) -> str:
        """The original tables this coverage was drawn from, rendered faithfully."""
        blocks: list[str] = []
        for result in results if results is not None else getattr(self, "_results", []):
            blocks.append(f"<h4>{escape(result.chapter_id)} :: {escape(result.caption)}</h4>")
            blocks.append(result.to_html())
        return "\n".join(blocks)

    def to_html(self, *, include_source: bool = False) -> str:
        """Render the gene's testing-coverage table (summary + detail)."""
        if not self.rows:
            missing = ", ".join(e.chapter_id for e in self.entries_indexed_only) or "none"
            return (
                f'<p class="gr-coverage-empty">No testing-strategy table found for '
                f"<strong>{escape(self.gene)}</strong>. "
                f"Indexed in {len(self.entries)} entr(ies); indexed without a table: "
                f"{escape(missing)}.</p>"
            )

        heading = f"{self.gene}"
        if self.gene_ids:
            heading += f" (NCBI Gene ID {', '.join(self.gene_ids)})"

        parts = [
            '<div class="gr-coverage">',
            f'<h3 class="gr-coverage-title">Testing coverage: {escape(heading)}</h3>',
        ]

        # Summary: one row per method.
        parts.append('<table class="gr-table gr-coverage-summary">')
        parts.append(
            "<caption>Expected yield by method</caption>"
            "<thead><tr>"
            '<th scope="col">Method</th><th scope="col">Kind</th>'
            '<th scope="col">Entries</th><th scope="col">Yield (observed)</th>'
            '<th scope="col">Range</th>'
            "</tr></thead><tbody>"
        )
        for row in self.summary_rows():
            if row["min_percent"] is None:
                rng = "&mdash;"
            elif row["min_percent"] == row["max_percent"]:
                rng = f"{row['min_percent']:g}%"
            else:
                rng = f"{row['min_percent']:g}&ndash;{row['max_percent']:g}%"
            if row["qualitative"]:
                rng += (
                    f'<br><span class="gr-qualitative">qualitative: '
                    f'{escape("; ".join(row["qualitative"]))}</span>'
                )
            parts.append(
                "<tr>"
                f'<th scope="row">{escape(row["method"])}</th>'
                f'<td>{escape(row["method_kind"])}</td>'
                f'<td>{row["entries"]}</td>'
                f'<td>{escape("; ".join(row["yields"]))}</td>'
                f"<td>{rng}</td>"
                "</tr>"
            )
        parts.append("</tbody></table>")

        # Detail: one row per entry x method.
        parts.append('<table class="gr-table gr-coverage-detail">')
        parts.append(
            "<caption>Per-entry detail</caption>"
            "<thead><tr>"
            '<th scope="col">Entry</th><th scope="col">Gene as written</th>'
            '<th scope="col">Method</th><th scope="col">Yield</th>'
            '<th scope="col">Status</th><th scope="col">Notes</th>'
            "</tr></thead><tbody>"
        )
        for row in self.rows:
            notes = []
            if row.order == 0:
                notes.append("listed first")
            if row.implicit:
                notes.append("gene inferred from chapter")
            if row.attribution:
                notes.append(f"gene attribution {row.attribution}")
            if row.detail:
                notes.append(row.detail)
            parts.append(
                "<tr>"
                f'<th scope="row">{escape(row.chapter_id)}</th>'
                f"<td>{escape(row.gene_label)}</td>"
                f"<td>{escape(row.method)}</td>"
                f"<td>{escape(row.yield_raw)}</td>"
                f"<td>{escape(row.observation.kind)}</td>"
                f'<td>{escape("; ".join(notes))}</td>'
                "</tr>"
            )
        parts.append("</tbody></table>")

        gaps = self.entries_indexed_only
        if gaps:
            names = ", ".join(f"{e.chapter_id}" for e in gaps)
            parts.append(
                f'<p class="gr-coverage-gap">Also indexed in {len(gaps)} entr(ies) with no '
                f"testing table: {escape(names)}.</p>"
            )

        if include_source:
            parts.append(self.source_html())
            parts.append("</div>")
        else:
            parts.append("</div>")
        return "\n".join(parts)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"<GeneTestingCoverage {self.gene!r} entries={len(self.entries)} "
            f"rows={len(self.rows)}>"
        )


# ──────────────────────────────────────────────────────────────────────────────
# Extraction
# ──────────────────────────────────────────────────────────────────────────────
def _table_has_gene_column(result: TestingYieldTable) -> bool:
    return any(cell_genes(row.gene) for row in result.rows)


def gene_testing_coverages(
    genes: Iterable[str],
    chapters: Iterable[Chapter],
    *,
    min_confidence: float = 0.8,
    require_caption_match: bool = False,
    include_indexed_without_table: bool = True,
) -> dict[str, GeneTestingCoverage]:
    """Build testing-coverage tables for several genes in a **single** corpus scan.

    Scanning the corpus dominates the cost (~5 s for 996 chapters), so querying a panel of
    genes through this function is dramatically cheaper than calling
    :func:`gene_testing_coverage` once per gene.

    Parameters
    ----------
    genes:
        Gene symbols, matched case-insensitively against the symbols inside a gene cell
        (so ``HBA2`` matches a row labelled ``HBA1 & HBA2``).
    chapters:
        Chapters to scan, e.g. ``corpus.iter_chapters()``.
    min_confidence:
        Drop weak table matches before looking for the genes (default 0.8).
    require_caption_match:
        Restrict to tables whose caption names genetic testing.
    include_indexed_without_table:
        Record entries that index a gene but publish no method/yield table, which is how
        coverage gaps stay visible (default ``True``).
    """
    targets = [g.strip().upper() for g in genes if g and g.strip()]
    if not targets:
        raise ValueError("at least one non-empty gene symbol is required")

    coverages = {t: GeneTestingCoverage(gene=t) for t in targets}
    sources: dict[str, list[TestingYieldTable]] = {t: [] for t in targets}

    for chapter in chapters:
        indexed = {g.upper() for g in chapter.genes}
        results = extract_testing_yield(
            chapter,
            require_caption_match=require_caption_match,
            min_confidence=min_confidence,
        )
        id_map = chapter.gene_id_map
        entries = {
            t: CoverageEntry(
                chapter_id=chapter.id, title=chapter.title, indexed=t in indexed
            )
            for t in targets
        }
        for t in targets:
            gene_id = id_map.get(t)
            if gene_id:
                entries[t].gene_id = gene_id
                if gene_id not in coverages[t].gene_ids:
                    coverages[t].gene_ids.append(gene_id)

        if results:
            for result in results:
                # Implicit attribution: a testing table with no usable Gene column belongs
                # to the chapter's single indexed gene, if it has exactly one.
                has_gene_column = _table_has_gene_column(result)
                for t in targets:
                    implicit = (not has_gene_column) and len(indexed) == 1 and t in indexed
                    if not has_gene_column and not implicit:
                        continue
                    matched = 0
                    for order, row in enumerate(result.rows):
                        if not implicit and t not in cell_genes(row.gene):
                            continue
                        coverages[t].rows.append(
                            CoverageRow(
                                gene=t,
                                gene_label=row.gene or t,
                                chapter_id=chapter.id,
                                chapter_title=chapter.title,
                                table_id=result.table_id,
                                caption=result.caption,
                                method=row.method,
                                method_raw=row.method_raw,
                                method_kind=row.method_kind,
                                is_sequencing=row.is_sequencing,
                                is_copy_number=row.is_copy_number,
                                yield_raw=row.yield_raw,
                                observation=row.observation,
                                attribution=row.attribution,
                                detail=row.detail,
                                order=order,
                                implicit=implicit,
                            )
                        )
                        matched += 1
                    entry = entries[t]
                    if matched:
                        entry.row_count += matched
                        if not entry.has_table:
                            entry.has_table = True
                            entry.table_id = result.table_id
                            entry.layout = result.layout
                            entry.confidence = result.confidence
                            sources[t].append(result)

        for t in targets:
            entry = entries[t]
            if entry.row_count or (entry.indexed and include_indexed_without_table):
                coverages[t].entries.append(entry)

    for t in targets:
        coverage = coverages[t]
        # Entries carrying data come first; then indexed-but-no-table gaps.
        coverage.entries.sort(key=lambda e: (not e.has_table, not e.indexed, e.chapter_id))
        coverage._results = sources[t]  # type: ignore[attr-defined]
    return coverages


def gene_testing_coverage(
    gene: str, chapters: Iterable[Chapter], **kwargs: object
) -> GeneTestingCoverage:
    """Build the testing-coverage table for a single ``gene``.

    See :func:`gene_testing_coverages` for the parameters. Use that function instead when
    querying several genes, so the corpus is scanned only once.
    """
    target = gene.strip().upper()
    if not target:
        raise ValueError("gene must be a non-empty symbol")
    return gene_testing_coverages([target], chapters, **kwargs)[target]


def coverage_matrix_html(coverages: Iterable[GeneTestingCoverage]) -> str:
    """Render several genes as one method x gene comparison table."""
    items = list(coverages)
    methods: list[str] = []
    for coverage in items:
        for row in coverage.summary_rows():
            if row["method"] not in methods:
                methods.append(row["method"])

    parts = [
        '<table class="gr-table gr-coverage-matrix">',
        "<caption>Testing coverage by gene</caption>",
        '<thead><tr><th scope="col">Gene</th>',
    ]
    parts.extend(f'<th scope="col">{escape(m)}</th>' for m in methods)
    parts.append("</tr></thead><tbody>")

    for coverage in items:
        parts.append(f'<tr><th scope="row">{escape(coverage.gene)}</th>')
        by_method = {r["method"]: r for r in coverage.summary_rows()}
        for method in methods:
            row = by_method.get(method)
            parts.append(
                "<td>&mdash;</td>" if not row else f'<td>{escape("; ".join(row["yields"]))}</td>'
            )
        parts.append("</tr>")
    parts.append("</tbody></table>")
    return "\n".join(parts)
