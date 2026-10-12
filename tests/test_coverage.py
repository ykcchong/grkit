"""Tests for :mod:`grkit.coverage` -- gene -> testing-coverage tables."""

from __future__ import annotations

import pytest

from grkit.coverage import (
    cell_genes,
    coverage_matrix_html,
    gene_testing_coverage,
    gene_testing_coverages,
)
from grkit.reader import Chapter


@pytest.fixture(scope="module")
def coverage_chapter(fixture_dir):
    return Chapter.from_file(str(fixture_dir / "coverage.nxml"))


@pytest.fixture(scope="module")
def solo_chapter(fixture_dir):
    return Chapter.from_file(str(fixture_dir / "single-gene.nxml"))


# ──────────────────────────────────────────────────────────────────────────────
# Gene-cell parsing
# ──────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize(
    "label,expected",
    [
        ("TEST1", ["TEST1"]),
        ("TEST1 & TEST2", ["TEST1", "TEST2"]),
        ("HBA1 & HBA2", ["HBA1", "HBA2"]),
        ("MT-TC MT-TF  MT-TH", ["MT-TC", "MT-TF", "MT-TH"]),
        ("CACNA1A, ATP1A2; SCN1A", ["CACNA1A", "ATP1A2", "SCN1A"]),
        ("PRS/AARS", ["PRS", "AARS"]),
        ("Unknown", []),
        ("NA", []),
        ("None reported", []),
        ("-", []),
        ("", []),
        ("SCN1A-AS1", ["SCN1A-AS1"]),
        # Symbols whose lowercase "orf" is legitimate.
        ("C11orf65", ["C11ORF65"]),
        ("C2orf71", ["C2ORF71"]),
        # Sentence-like cells in a "Gene" column (recurrent-deletion entries) yield nothing.
        ("~0.8-Mb heterozygous deletion at 1q21.1", []),
        ("recurrent 1q21.1 to 1q21.2 hybrid allele", []),
        ("ISCN: 46,XY", []),
        ("SEQ[GRCH38]", []),
        ("ClinGen ID", []),
        ("ISCA-37392", []),
        # Assay abbreviations and chromosome arms are not genes.
        ("MLPA", []),
        ("1q21.1", []),
    ],
)
def test_cell_genes(label, expected):
    assert cell_genes(label) == expected


def test_cell_genes_keeps_a_long_mitochondrial_run():
    label = "MT-TC MT-TF MT-TH MT-TK MT-TL2 MT-TQ MT-TV MT-TW MT-TS1 MT-TS2 MT-ND1 MT-ND6"
    genes = cell_genes(label)
    assert len(genes) == 12
    assert genes[0] == "MT-TC" and genes[-1] == "MT-ND6"


def test_cell_genes_does_not_split_hyphenated_symbols():
    """A query for SCN1A must not match the antisense gene SCN1A-AS1."""
    assert "SCN1A" not in cell_genes("SCN1A-AS1")


# ──────────────────────────────────────────────────────────────────────────────
# Coverage extraction
# ──────────────────────────────────────────────────────────────────────────────
def test_multi_gene_cell_matches_either_gene(coverage_chapter):
    """Both TEST1 and TEST2 must find the shared "TEST1 & TEST2" row."""
    for gene in ("TEST1", "TEST2"):
        coverage = gene_testing_coverage(gene, [coverage_chapter])
        assert coverage.has_coverage
        labels = {row.gene_label for row in coverage.rows}
        assert any("&" in label for label in labels), f"{gene} missed the shared row"


def test_coverage_rows_for_test1(coverage_chapter):
    coverage = gene_testing_coverage("TEST1", [coverage_chapter])
    observed = [(r.gene_label, r.method, r.yield_raw) for r in coverage.rows]
    assert observed == [
        ("TEST1 & TEST2", "Sequence analysis", "50%"),
        ("TEST1", "Gene-targeted deletion/duplication analysis", "25%"),
    ]
    assert all(r.chapter_id == "cov" for r in coverage.rows)
    assert all(r.implicit is False for r in coverage.rows)


def test_gene_ids_come_from_the_genes_and_databases_table(coverage_chapter):
    assert gene_testing_coverage("TEST1", [coverage_chapter]).gene_ids == ["111"]
    assert gene_testing_coverage("TEST2", [coverage_chapter]).gene_ids == ["222"]


def test_gene_id_map_only_covers_published_genes(coverage_chapter):
    assert coverage_chapter.gene_id_map == {"TEST1": "111", "TEST2": "222"}
    assert coverage_chapter.gene_ids == ["111", "222"]


def test_coverage_gap_is_reported_not_dropped(coverage_chapter):
    """ABSENT is indexed but publishes no testing table."""
    coverage = gene_testing_coverage("ABSENT", [coverage_chapter])
    assert coverage.has_coverage is False
    assert coverage.rows == []
    assert [e.chapter_id for e in coverage.entries_indexed_only] == ["cov"]
    assert coverage.entries[0].indexed is True
    assert coverage.entries[0].has_table is False


def test_unknown_gene_yields_empty_coverage(coverage_chapter):
    coverage = gene_testing_coverage("NOSUCHGENE", [coverage_chapter])
    assert coverage.entries == []
    assert coverage.rows == []


def test_empty_gene_is_rejected(coverage_chapter):
    with pytest.raises(ValueError):
        gene_testing_coverage("   ", [coverage_chapter])


# ──────────────────────────────────────────────────────────────────────────────
# Implicit attribution (table with no Gene column)
# ──────────────────────────────────────────────────────────────────────────────
def test_single_gene_chapter_attributes_rows_implicitly(solo_chapter):
    coverage = gene_testing_coverage("SOLO", [solo_chapter])
    assert [r.method for r in coverage.rows] == [
        "Sequence analysis",
        "Gene-targeted deletion/duplication analysis",
    ]
    assert [r.yield_raw for r in coverage.rows] == ["88%", "12%"]
    assert all(r.implicit is True for r in coverage.rows)
    assert coverage.gene_ids == ["999"]


def test_implicit_attribution_does_not_leak_to_other_genes(solo_chapter):
    """A gene the chapter does not index must not inherit its rows."""
    coverage = gene_testing_coverage("OTHER", [solo_chapter])
    assert coverage.rows == []


def test_implicit_rows_are_flagged_in_html(solo_chapter):
    html = gene_testing_coverage("SOLO", [solo_chapter]).to_html()
    assert "gene inferred from chapter" in html


# ──────────────────────────────────────────────────────────────────────────────
# Aggregation
# ──────────────────────────────────────────────────────────────────────────────
def test_summary_rows(coverage_chapter):
    coverage = gene_testing_coverage("TEST1", [coverage_chapter])
    summary = {row["method"]: row for row in coverage.summary_rows()}
    sequence = summary["Sequence analysis"]
    assert sequence["method_kind"] == "sequence"
    assert sequence["is_sequencing"] is True
    assert sequence["entries"] == 1
    assert sequence["yields"] == ["50%"]
    assert sequence["min_percent"] == sequence["max_percent"] == 50.0
    assert sequence["first_order"] == 0


def test_qualitative_yields_are_kept_out_of_the_numeric_range(coverage_chapter):
    coverage = gene_testing_coverage("TEST2", [coverage_chapter])
    summary = {
        row["method"]: row for row in coverage.summary_rows()
    }["Gene-targeted deletion/duplication analysis"]
    assert summary["yields"] == ["None reported"]
    # "None reported" is an explicit zero, so the range is 0 rather than unknown.
    assert summary["min_percent"] == 0.0
    assert summary["qualitative"] == []


def test_sensitivity_and_preferred_strategy(coverage_chapter):
    coverage = gene_testing_coverage("TEST1", [coverage_chapter])
    sensitivity = coverage.sensitivity()
    assert sensitivity["sequencing"]["yields"] == ["50%"]
    assert sensitivity["copy_number"]["yields"] == ["25%"]
    assert sensitivity["sequencing"]["best_percent"] == 50.0
    assert sensitivity["preferred"] == "sequence", "GeneReviews lists the first-line method first"


def test_first_line_is_the_first_listed_method(coverage_chapter):
    first = gene_testing_coverage("TEST1", [coverage_chapter]).first_line()
    assert first is not None
    assert first.method == "Sequence analysis"
    assert first.order == 0


def test_coverage_matrix(coverage_chapter):
    matrix = gene_testing_coverage("TEST1", [coverage_chapter]).coverage_matrix()
    assert matrix == {
        "cov": {
            "Sequence analysis": "50%",
            "Gene-targeted deletion/duplication analysis": "25%",
        }
    }


# ──────────────────────────────────────────────────────────────────────────────
# Multi-gene scan
# ──────────────────────────────────────────────────────────────────────────────
def test_multi_gene_scan_matches_individual_queries(coverage_chapter):
    many = gene_testing_coverages(["TEST1", "TEST2", "ABSENT"], [coverage_chapter])
    assert set(many) == {"TEST1", "TEST2", "ABSENT"}
    for gene in many:
        single = gene_testing_coverage(gene, [coverage_chapter])
        assert [r.to_dict() for r in many[gene].rows] == [r.to_dict() for r in single.rows]
        assert [e.to_dict() for e in many[gene].entries] == [
            e.to_dict() for e in single.entries
        ]


def test_multi_gene_scan_normalises_case(coverage_chapter):
    many = gene_testing_coverages(["test1"], [coverage_chapter])
    assert many["TEST1"].gene == "TEST1"


def test_entries_with_data_sort_first(coverage_chapter):
    coverage = gene_testing_coverage("ABSENT", [coverage_chapter])
    assert coverage.entries  # gap entry present
    combined = gene_testing_coverage("TEST1", [coverage_chapter])
    assert combined.entries_with_table


# ──────────────────────────────────────────────────────────────────────────────
# Serialization and rendering
# ──────────────────────────────────────────────────────────────────────────────
def test_to_dict_shape(coverage_chapter):
    record = gene_testing_coverage("TEST1", [coverage_chapter]).to_dict()
    assert record["gene"] == "TEST1"
    assert record["gene_ids"] == ["111"]
    assert record["entries_with_table"] == 1
    assert record["has_coverage"] is True
    assert record["rows"][0]["yield_percent"] == 50.0
    assert record["matrix"]["cov"]["Sequence analysis"] == "50%"


def test_to_records(coverage_chapter):
    records = gene_testing_coverage("TEST1", [coverage_chapter]).to_records()
    assert len(records) == 2
    assert {"gene", "chapter_id", "method", "yield_raw", "implicit"} <= set(records[0])


def test_to_html_contains_summary_and_detail(coverage_chapter):
    html = gene_testing_coverage("TEST1", [coverage_chapter]).to_html()
    assert "Testing coverage: TEST1 (NCBI Gene ID 111)" in html
    assert "Expected yield by method" in html
    assert "Per-entry detail" in html
    assert "Sequence analysis" in html
    assert "listed first" in html
    assert "gr-coverage-summary" in html


def test_to_html_reports_gaps(coverage_chapter):
    html = gene_testing_coverage("ABSENT", [coverage_chapter]).to_html()
    assert "No testing-strategy table found for" in html
    assert "cov" in html


def test_to_html_appends_source_tables(coverage_chapter):
    html = gene_testing_coverage(
        "TEST1", [coverage_chapter]
    ).to_html(include_source=True)
    assert "<table" in html
    assert "Test Coverage Disorder" in html


def test_coverage_matrix_html(coverage_chapter):
    coverages = gene_testing_coverages(["TEST1", "TEST2"], [coverage_chapter])
    html = coverage_matrix_html(coverages.values())
    assert "Testing coverage by gene" in html
    assert "gr-coverage-matrix" in html
    assert html.count("<tr>") == 3  # header + two genes


# ──────────────────────────────────────────────────────────────────────────────
# Corpus integration
# ──────────────────────────────────────────────────────────────────────────────
def test_corpus_gene_coverage(tmp_path, fixture_dir):
    import shutil

    from grkit import Corpus

    shutil.copy(fixture_dir / "coverage.nxml", tmp_path / "cov.nxml")
    with Corpus(tmp_path) as corpus:
        coverage = corpus.gene_coverage("TEST2")
        assert coverage.gene == "TEST2"
        assert len(coverage.rows) == 2
        many = corpus.gene_coverages(["TEST1", "TEST2"])
        assert set(many) == {"TEST1", "TEST2"}
