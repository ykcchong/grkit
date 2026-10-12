"""Tests for :mod:`grkit.sources` and :mod:`grkit.corpus`."""

from __future__ import annotations

import io
import os
import tarfile

import pytest

from grkit.corpus import Corpus
from grkit.sources import (
    DirectorySource,
    FileSource,
    TarGzSource,
    is_archive,
    open_source,
)

CHAPTER = """<?xml version="1.0" encoding="UTF-8"?>
<book-part-wrapper id="{cid}" content-type="chapter">
  <book-part><book-part-meta>
    <title-group><title>{title}</title></title-group>
    <kwd-group kwd-group-type="GeneSymbol"><kwd>{gene}</kwd></kwd-group>
    <kwd-group kwd-group-type="DiseaseName"><kwd>{title}</kwd></kwd-group>
  </book-part-meta>
  <body>
    <sec id="{cid}.Diagnosis"><title>Diagnosis</title><p>Body text for {title}.</p></sec>
  </body></book-part>
</book-part-wrapper>
"""

MALFORMED = "<book-part-wrapper><this is not xml"


def make_chapter(cid: str, title: str = "Test Disorder", gene: str = "TEST1") -> str:
    return CHAPTER.format(cid=cid, title=title, gene=gene)


@pytest.fixture
def corpus_dir(tmp_path):
    root = tmp_path / "nxmls"
    root.mkdir()
    (root / "alpha.nxml").write_text(make_chapter("alpha", "Alpha Syndrome", "AAA"))
    (root / "beta.nxml").write_text(make_chapter("beta", "Beta Disorder", "BBB"))
    (root / "broken.nxml").write_text(MALFORMED)
    (root / "notes.txt").write_text("not a chapter")
    return root


def build_targz(path, members: dict[str, str]) -> str:
    with tarfile.open(path, "w:gz") as tar:
        for name, content in members.items():
            data = content.encode("utf-8")
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return str(path)


@pytest.fixture
def archive(tmp_path):
    return build_targz(
        tmp_path / "gene_NBK1116.tar.gz",
        {
            "gene_NBK1116/alpha.nxml": make_chapter("alpha", "Alpha Syndrome", "AAA"),
            "gene_NBK1116/beta.nxml": make_chapter("beta", "Beta Disorder", "BBB"),
            "gene_NBK1116/ibm-Image001.jpg": "binary-ish",
            "gene_NBK1116/README.txt": "ignore me",
        },
    )


# ──────────────────────────────────────────────────────────────────────────────
# DirectorySource / FileSource
# ──────────────────────────────────────────────────────────────────────────────
def test_directory_source_lists_only_nxml(corpus_dir):
    source = DirectorySource(corpus_dir)
    assert source.names() == ["alpha.nxml", "beta.nxml", "broken.nxml"]
    assert len(source) == 3
    assert "alpha.nxml" in source and "notes.txt" not in source


def test_directory_source_reads_with_or_without_suffix(corpus_dir):
    source = DirectorySource(corpus_dir)
    assert source.read("alpha.nxml") == source.read("alpha")
    assert b"book-part-wrapper" in source.read("alpha")


def test_directory_source_iterates_in_order(corpus_dir):
    names = [name for name, _ in DirectorySource(corpus_dir).iter_bytes()]
    assert names == ["alpha.nxml", "beta.nxml", "broken.nxml"]


def test_file_source(corpus_dir):
    source = FileSource(corpus_dir / "alpha.nxml")
    assert source.names() == ["alpha.nxml"]
    assert source.read("alpha") == source.read("alpha.nxml")
    with pytest.raises(KeyError):
        source.read("beta")


def test_directory_source_rejects_a_file(tmp_path):
    path = tmp_path / "x.nxml"
    path.write_text(make_chapter("x"))
    with pytest.raises(NotADirectoryError):
        DirectorySource(path)


# ──────────────────────────────────────────────────────────────────────────────
# TarGzSource
# ──────────────────────────────────────────────────────────────────────────────
def test_targz_names_skip_non_chapters(archive):
    assert TarGzSource(archive).names() == ["alpha.nxml", "beta.nxml"]


def test_targz_streaming_yields_every_chapter(archive):
    with TarGzSource(archive) as source:
        names = [name for name, _ in source.iter_bytes()]
    assert names == ["alpha.nxml", "beta.nxml"]


def test_targz_random_access_by_basename(archive):
    with TarGzSource(archive) as source:
        assert b"Alpha Syndrome" in source.read("alpha.nxml")
        assert source.read("beta") == source.read("beta.nxml")


def test_targz_missing_member_raises_keyerror(archive):
    with TarGzSource(archive) as source:
        with pytest.raises(KeyError):
            source.read("nope")


def test_targz_scan_budget_defers_spooling(archive):
    """A single lookup should not materialise the whole tar."""
    source = TarGzSource(archive, scan_budget=3)
    try:
        assert b"Alpha Syndrome" in source.read("alpha")
        assert not source.spooled
        assert source._spool_path is None
    finally:
        source.close()


def test_targz_spools_after_budget_and_cleans_up(archive):
    source = TarGzSource(archive, scan_budget=0)
    try:
        assert b"Alpha Syndrome" in source.read("alpha")
        assert source.spooled
        spool_path = source._spool_path
        assert spool_path and os.path.exists(spool_path)
    finally:
        source.close()
    assert not os.path.exists(spool_path), "the temporary tar must be deleted on close"


def test_targz_spool_disabled_writes_nothing(archive):
    source = TarGzSource(archive, spool=False)
    try:
        assert b"Alpha Syndrome" in source.read("alpha")
        assert not source.spooled
        assert source._spool_path is None
    finally:
        source.close()


def test_targz_handles_uncompressed_tar(tmp_path):
    plain = tmp_path / "corpus.tar"
    with tarfile.open(plain, "w") as tar:
        data = make_chapter("solo").encode()
        info = tarfile.TarInfo("solo.nxml")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    with TarGzSource(plain) as source:
        assert source.names() == ["solo.nxml"]
        assert b"solo" in source.read("solo")


def test_targz_accepts_a_file_object(archive):
    """The stream is also readable through the generic ``open`` helper."""
    with TarGzSource(archive) as source:
        with source.open("alpha.nxml") as handle:
            assert b"Alpha Syndrome" in handle.read()


# ──────────────────────────────────────────────────────────────────────────────
# Dispatch
# ──────────────────────────────────────────────────────────────────────────────
def test_is_archive_detects_forms(corpus_dir, archive):
    assert is_archive(archive)
    assert not is_archive(corpus_dir)
    assert not is_archive(corpus_dir / "alpha.nxml")


@pytest.mark.parametrize("kind", ["dir", "archive", "file"])
def test_open_source_dispatch(kind, corpus_dir, archive):
    target = {"dir": corpus_dir, "archive": archive, "file": corpus_dir / "alpha.nxml"}[kind]
    expected = {"dir": DirectorySource, "archive": TarGzSource, "file": FileSource}[kind]
    source = open_source(target)
    try:
        assert isinstance(source, expected)
        assert source.names()
    finally:
        source.close()


def test_open_source_missing_path(tmp_path):
    with pytest.raises(FileNotFoundError):
        open_source(tmp_path / "nope.tar.gz")


# ──────────────────────────────────────────────────────────────────────────────
# Corpus
# ──────────────────────────────────────────────────────────────────────────────
def test_corpus_over_directory(corpus_dir):
    corpus = Corpus(corpus_dir)
    assert len(corpus) == 3
    assert corpus.ids == ["alpha", "beta", "broken"]
    assert "alpha" in corpus and "nope" not in corpus


def test_corpus_get_and_missing(corpus_dir):
    corpus = Corpus(corpus_dir)
    chapter = corpus.get("alpha")
    assert chapter is not None and chapter.title == "Alpha Syndrome"
    assert corpus.get("alpha.nxml") is not None
    assert corpus.get("does-not-exist") is None
    assert corpus.get("broken") is None, "malformed XML yields None, not an exception"


def test_corpus_iter_skips_malformed_by_default(corpus_dir):
    corpus = Corpus(corpus_dir)
    assert [c.id for c in corpus] == ["alpha", "beta"]


def test_corpus_iter_can_raise(corpus_dir):
    import xml.etree.ElementTree as ET

    corpus = Corpus(corpus_dir)
    with pytest.raises(ET.ParseError):
        list(corpus.iter_chapters(skip_errors=False))


def test_corpus_find_by_gene_and_disease(corpus_dir):
    corpus = Corpus(corpus_dir)
    assert [c.id for c in corpus.find_by_gene("aaa")] == ["alpha"]
    assert [c.id for c in corpus.find_by_gene("BBB")] == ["beta"]
    assert [c.id for c in corpus.find_by_disease("Beta")] == ["beta"]
    assert corpus.find_by_gene("zzz") == []


def test_corpus_over_archive_streams(archive):
    corpus = Corpus(archive)
    assert len(corpus) == 2
    assert [c.id for c in corpus] == ["alpha", "beta"]
    assert corpus.get("beta").genes == ["BBB"]


def test_corpus_search_matches_body_text(corpus_dir):
    corpus = Corpus(corpus_dir)
    assert [c.id for c in corpus.search("Body text for Alpha")] == ["alpha"]


def test_corpus_iter_tables_and_testing(corpus_dir):
    corpus = Corpus(corpus_dir)
    assert list(corpus.iter_tables()) == []
    assert list(corpus.iter_testing_yield()) == []


def test_corpus_closes_only_owned_sources(corpus_dir):
    source = DirectorySource(corpus_dir)
    Corpus(source).close()
    assert source.names(), "an externally supplied source must not be closed"
