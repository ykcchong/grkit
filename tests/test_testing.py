"""Tests for :mod:`grkit.testing` -- the testing-strategy-vs-yield extractor."""

from __future__ import annotations

import csv
import json

import pytest

from grkit.testing import (
    classify_method,
    extract_from_table,
    extract_testing_yield,
    is_testing_yield_table,
    parse_yield,
    to_records,
    write_csv,
    write_json,
)


@pytest.fixture
def results(chapter):
    return {r.table_id: r for r in extract_testing_yield(chapter)}


# ──────────────────────────────────────────────────────────────────────────────
# Table selection
# ──────────────────────────────────────────────────────────────────────────────
def test_only_testing_tables_are_selected(results):
    assert set(results) == {"tst.T.testing", "tst.T.wide", "tst.T.twoyield"}


def test_gene_database_table_is_not_a_testing_table(table_of):
    """'Chromosome Locus' is a data descriptor, not a cytogenetic test method."""
    assert not is_testing_yield_table(table_of("tst.T.db"))
    assert extract_from_table(table_of("tst.T.db")) is None


def test_confidence_and_caption_match(results):
    for result in results.values():
        assert result.caption_match is True
        assert result.confidence >= 0.9
        assert result.title == "Test Syndrome"


# ──────────────────────────────────────────────────────────────────────────────
# Long layout
# ──────────────────────────────────────────────────────────────────────────────
def test_long_layout_pairs_method_with_yield(results):
    result = results["tst.T.testing"]
    assert result.layout == "long"
    observed = [
        (r.gene, r.method, r.yield_raw, r.method_kind) for r in result.rows
    ]
    assert observed == [
        ("TEST1", "Sequence analysis", "92%", "sequence"),
        ("TEST1", "Gene-targeted deletion/duplication analysis", "2.5%", "deletion_duplication"),
        ("Unknown", "NA", "5.5%", "other"),
    ]


def test_long_layout_flags_sequencing_and_copy_number(results):
    rows = results["tst.T.testing"].rows
    assert [r.is_sequencing for r in rows] == [True, False, False]
    assert [r.is_copy_number for r in rows] == [False, True, False]


def test_strategy_vs_yield_comparison(results):
    """The headline view: sequencing and copy-number yield side by side."""
    comparison = {c["gene"]: c for c in results["tst.T.testing"].comparison()}
    assert comparison["TEST1"]["sequencing"] == "92%"
    assert comparison["TEST1"]["sequencing_percent"] == 92.0
    assert comparison["TEST1"]["copy_number"] == "2.5%"
    assert comparison["TEST1"]["copy_number_percent"] == 2.5
    assert comparison["Unknown"]["other_methods"] == [
        {
            "method": "NA",
            "yield": "5.5%",
            "column": "Proportion of Probands with a Pathogenic Variant Detectable by Method",
        }
    ]


def test_footnote_marks_are_captured_on_rows(results):
    row = next(
        r for r in results["tst.T.testing"].rows if r.method_kind == "deletion_duplication"
    )
    assert set(row.marks) == {"4", "5"}
    assert row.yield_raw == "2.5%"


# ──────────────────────────────────────────────────────────────────────────────
# Wide layout (methods as columns)
# ──────────────────────────────────────────────────────────────────────────────
def test_wide_layout_methods_come_from_the_header_path(results):
    result = results["tst.T.wide"]
    assert result.layout == "wide"
    assert result.methods == [
        "Sequence analysis",
        "Gene-targeted deletion/duplication analysis – Common deletions",
        "Gene-targeted deletion/duplication analysis – Other deletions",
    ]


def test_wide_layout_yields_and_attribution(results):
    result = results["tst.T.wide"]
    rows = [r for r in result.rows if r.gene == "TEST1"]
    assert [r.yield_raw for r in rows] == ["~15%", "~85%", "<5%"]
    assert all(r.attribution == ">98%" for r in rows)
    assert all(r.is_copy_number for r in rows[1:])


def test_wide_layout_nested_submethods_stay_copy_number(results):
    rows = [r for r in results["tst.T.wide"].rows if r.method_kind == "deletion_duplication"]
    assert rows, "nested deletion/duplication sub-methods must classify as copy-number"


def test_two_yield_columns_are_distinguished(results):
    result = results["tst.T.twoyield"]
    assert {r.yield_column for r in result.rows} == {"Proband", "At-risk family members"}
    assert len(result.rows) == 4
    assert {r.method for r in result.rows} == {"CMA", "Targeted deletion analysis"}


def test_comparison_prefers_a_populated_submethod(results):
    """A blank nested cell must not shadow a sibling that states a yield."""
    comparison = {c["gene"]: c for c in results["tst.T.twoyield"].comparison()}
    # The gene key is the 'Deletion' column value.
    gene = next(iter(comparison))
    assert comparison[gene]["copy_number"] == "100%"


# ──────────────────────────────────────────────────────────────────────────────
# Serialization
# ──────────────────────────────────────────────────────────────────────────────
def test_to_dict_and_source_html(results):
    record = results["tst.T.testing"].to_dict(include_source_html=True)
    assert record["layout"] == "long"
    assert record["has_sequencing"] and record["has_copy_number"]
    assert record["strategy_vs_yield"][0]["gene"] == "TEST1"
    assert 'rowspan="2"' in record["source_table_html"]


def test_flat_records_include_table_context(results):
    records = to_records(results.values())
    assert records
    first = records[0]
    assert {"chapter_id", "table_id", "caption", "gene", "method", "yield_raw"} <= set(first)


def test_write_json_and_csv(results, tmp_path):
    payload = tmp_path / "out.json"
    write_json(results.values(), str(payload))
    data = json.loads(payload.read_text())
    assert len(data) == 3
    assert data[0]["rows"]

    flat = tmp_path / "out.csv"
    write_csv(results.values(), str(flat))
    with open(flat, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == sum(len(r.rows) for r in results.values())


def test_require_caption_match_filters(chapter):
    assert len(extract_testing_yield(chapter, require_caption_match=True)) == 3


def test_min_confidence_filters(chapter):
    assert extract_testing_yield(chapter, min_confidence=0.99) == []


# ──────────────────────────────────────────────────────────────────────────────
# Yield parsing
# ──────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize(
    "raw,kind,point,low,high,qualifier",
    [
        ("92%", "percent", 92.0, 92.0, 92.0, ""),
        ("2.5%", "percent", 2.5, 2.5, 2.5, ""),
        ("0%", "percent", 0.0, 0.0, 0.0, ""),
        ("<1%", "percent", None, None, 1.0, "<"),
        (">99%", "percent", None, 99.0, None, ">"),
        ("~100%", "percent", 100.0, 100.0, 100.0, "~"),
        ("Nearly 100%", "percent", 100.0, 100.0, 100.0, "~"),
        ("Almost 100%", "percent", 100.0, 100.0, 100.0, "~"),
        ("Up to 100%", "percent", None, None, 100.0, "<="),
        ("Estimated >95%", "percent", None, 95.0, None, ">"),
        ("greater than 50%", "percent", None, 50.0, None, ">"),
        ("~70%-80%", "percent", None, 70.0, 80.0, "~"),
        ("5-15%", "percent", None, 5.0, 15.0, ""),
        ("15/15", "fraction", 100.0, 100.0, 100.0, ""),
        ("2/2 families", "fraction", 100.0, 100.0, 100.0, ""),
        ("129/139 (93%)", "percent", 93.0, 93.0, 93.0, ""),
        ("100% (14/14 families)", "percent", 100.0, 100.0, 100.0, ""),
        ("~97.5% (78/80)", "percent", 97.5, 97.5, 97.5, "~"),
        ("2.9% of sporadic cases 6.3% of familial cases", "multi_percent", None, 2.9, 2.9, ""),
        ("100", "percent", 100.0, 100.0, 100.0, ""),
    ],
)
def test_parse_yield_values(raw, kind, point, low, high, qualifier):
    observation = parse_yield(raw)
    assert observation.kind == kind
    assert observation.point == point
    assert observation.low == low
    assert observation.high == high
    assert observation.qualifier == qualifier


@pytest.mark.parametrize(
    "raw,kind",
    [
        ("NA", "not_applicable"),
        ("N/A", "not_applicable"),
        ("Not applicable", "not_applicable"),
        ("Unknown", "unknown"),
        ("Unknown, none reported", "unknown"),
        ("none reported", "none_reported"),
        ("See footnote 7.", "see_footnote"),
        ("See footnotes 2, 3.", "see_footnote"),
        ("1 reported", "count"),
        ("2 individuals", "count"),
        ("", "empty"),
        ("Rare", "text"),
        ("Normal", "text"),
    ],
)
def test_parse_yield_non_numeric(raw, kind):
    assert parse_yield(raw).kind == kind


def test_parse_yield_fraction_attaches_denominator():
    observation = parse_yield("100% (14/14 families)")
    assert observation.numerator == 14
    assert observation.denominator == 14


def test_parse_yield_counts_are_not_percentages():
    observation = parse_yield("3 families")
    assert observation.kind == "count"
    assert observation.count == 3.0
    assert observation.percent is None


def test_parse_yield_bare_number_as_count():
    observation = parse_yield("5", assume_percent=False)
    assert observation.kind == "count"
    assert observation.count == 5.0


def test_parse_yield_midpoint_uses_bound():
    assert parse_yield("<1%").midpoint == 1.0
    assert parse_yield(">98%").midpoint == 98.0
    assert parse_yield("~70%-80%").midpoint == 75.0


# ──────────────────────────────────────────────────────────────────────────────
# Method classification
# ──────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize(
    "label,kind,sequencing,copy_number",
    [
        ("Sequence analysis", "sequence", True, False),
        ("Sequence analysis of select exons", "sequence", True, False),
        ("Gene-targeted deletion/duplication analysis", "deletion_duplication", False, True),
        ("Targeted deletion analysis", "deletion_duplication", False, True),
        ("Targeted duplication analysis", "deletion_duplication", False, True),
        ("MLPA", "deletion_duplication", False, True),
        ("CMA", "cytogenetic", False, False),
        ("Chromosomal microarray analysis", "cytogenetic", False, False),
        ("Karyotype", "cytogenetic", False, False),
        ("FISH", "cytogenetic", False, False),
        ("Methylation analysis", "methylation", False, False),
        ("Repeat expansion testing", "repeat_expansion", False, False),
        ("Exome sequencing", "panel_or_genome", False, False),
        ("NA", "other", False, False),
        ("", "other", False, False),
    ],
)
def test_classify_method(label, kind, sequencing, copy_number):
    result = classify_method(label)
    assert result.kind == kind
    assert result.is_sequencing is sequencing
    assert result.is_copy_number is copy_number


def test_duplication_label_wins_over_targeted():
    """Order matters: del/dup must not be classified as generic 'targeted'."""
    assert classify_method("Gene-targeted deletion/duplication analysis").kind == (
        "deletion_duplication"
    )
