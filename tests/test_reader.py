"""Tests for :mod:`grkit.reader` -- chapter metadata, sections, references, records."""

from __future__ import annotations

import json

from grkit.reader import Chapter, classify_url, section_name


def test_identity_and_metadata(chapter):
    assert chapter.id == "tst"
    assert chapter.content_type == "chapter"
    assert chapter.title == "Test Syndrome"
    assert chapter.alt_titles == ["Synonyms: TS, Test Disorder"]
    assert chapter.genes == ["TEST1", "TEST2"]
    assert chapter.diseases == ["Test Syndrome"]
    assert chapter.synonyms == ["TS", "Test Disorder"]
    assert chapter.included_phenotypes == ["Test-related phenotype"]
    assert chapter.proteins == ["Test protein 1"]
    assert chapter.is_retired is False


def test_authors(chapter):
    authors = chapter.authors
    assert len(authors) == 1
    assert authors[0].full_name == "Jane Doe"
    assert authors[0].surname == "Doe"
    assert authors[0].email == "jane@example.org"


def test_dates(chapter):
    assert chapter.created_date == "2001-02-01"
    assert chapter.updated_date == "2020-04-03"


def test_sections_and_lookup(chapter):
    names = chapter.section_names
    assert names[0] == "Diagnosis"
    assert "Molecular_Genetics" in names
    assert chapter.get_section("Diagnosis") is not None
    assert chapter.get_section("diagnosis") is not None, "lookup is case-insensitive"
    assert chapter.get_section("tst.Diagnosis") is not None
    assert chapter.get_section("Nope") is None


def test_section_content_and_tables(chapter):
    diagnosis = chapter.get_section("Diagnosis")
    assert any("Table 1" in p for p in diagnosis.paragraphs)
    assert [t.id for t in diagnosis.tables] == ["tst.T.testing"]


def test_summary_is_structured(chapter):
    summary = chapter.summary
    assert summary[0]["title"] == "Clinical characteristics"
    assert "model disorder" in summary[0]["content"]


def test_references(chapter):
    refs = chapter.references
    assert len(refs) == 1
    ref = refs[0]
    assert ref.title == "A test article."
    assert ref.source == "J Test"
    assert ref.year == "2001"
    assert ref.pmid == "12345678"
    assert ref.authors == ["Doe J", "Roe R"]


def test_links_and_cross_references(chapter):
    hrefs = {link.href for link in chapter.links}
    assert "https://www.omim.org/entry/123456" in hrefs
    cross = {(x.db, x.value) for x in chapter.cross_references}
    assert ("omim", "123456") in cross


def test_cross_references_are_deduplicated(chapter):
    keys = [(x.db, x.value) for x in chapter.cross_references]
    assert len(keys) == len(set(keys))


def test_related_chapters(chapter):
    assert chapter.previous_chapter[0]["document_id"] == "prev-chap"
    assert chapter.next_chapter[0]["document_id"] == "next-chap"


def test_tables_are_tagged_with_owning_section(chapter):
    by_id = {t.id: t for t in chapter.tables}
    assert by_id["tst.T.testing"].section == "Diagnosis"
    assert by_id["tst.T.wide"].section == "Molecular_Genetics"
    assert all(t.chapter_id == "tst" for t in chapter.tables)


def test_table_by_id(chapter):
    assert chapter.table_by_id("tst.T.wide") is not None
    assert chapter.table_by_id("nope") is None


def test_section_tree_nesting(chapter):
    tree = chapter.section_tree()
    assert tree[0]["name"] == "Diagnosis"
    assert isinstance(tree[0]["subsections"], list)


def test_text_length_is_positive(chapter):
    assert chapter.text_length > 0


def test_to_record_is_json_serializable(chapter):
    record = chapter.to_record()
    text = json.dumps(record)  # must not raise
    assert '"title": "Test Syndrome"' in text
    assert record["genes"] == ["TEST1", "TEST2"]
    assert len(record["tables"]) == 5
    assert any("html" in table for table in record["tables"])


def test_to_record_can_omit_table_html(chapter):
    record = chapter.to_record(table_html=False)
    assert all("html" not in table for table in record["tables"])


def test_try_from_file_on_malformed_xml(tmp_path):
    bad = tmp_path / "bad.nxml"
    bad.write_text("<not xml")
    assert Chapter.try_from_file(str(bad)) is None


def test_from_string_roundtrip(sample_path):
    chapter = Chapter.from_string(sample_path.read_bytes())
    assert chapter.title == "Test Syndrome"
    assert chapter.id == "tst"


def test_retired_chapter_detected_from_title():
    xml = (
        '<book-part-wrapper id="old"><book-part><book-part-meta>'
        "<title-group><title>Old Disorder - RETIRED CHAPTER, FOR HISTORICAL "
        "REFERENCE ONLY</title></title-group></book-part-meta><body/></book-part>"
        "</book-part-wrapper>"
    )
    chapter = Chapter.from_string(xml)
    assert chapter.is_retired is True


def test_chapter_without_body_is_safe():
    xml = (
        '<book-part-wrapper id="x"><book-part><book-part-meta>'
        "<title-group><title>T</title></title-group></book-part-meta>"
        "</book-part></book-part-wrapper>"
    )
    chapter = Chapter.from_string(xml)
    assert chapter.sections == []
    assert chapter.tables == []
    assert chapter.text_length == 0


def test_section_name_helper():
    assert section_name("abs.Molecular_Genetics") == "Molecular_Genetics"
    assert section_name("Plain") == "Plain"


def test_classify_url_patterns():
    assert classify_url("https://www.omim.org/entry/123456") == ("omim", "123456")
    assert classify_url("https://www.ncbi.nlm.nih.gov/clinvar/variation/99/") == ("clinvar", "99")
    assert classify_url("https://pubmed.ncbi.nlm.nih.gov/12345/") == ("pmid", "12345")
    assert classify_url("https://example.org/whatever") == ("uri", "https://example.org/whatever")
    assert classify_url("") == ("", "")
