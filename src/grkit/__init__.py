"""grkit -- a robust reader for the GeneReviews (NCBI Bookshelf NBK1116) NXML corpus.

Designed around three requirements:

1. **Self-updating access.** :func:`~grkit.updater.fetch_archive` resolves the current
   archive URL from NCBI, downloads ``gene_NBK1116.tar.gz`` once, and caches it; the
   library reads chapters *directly* out of that tarball, so nothing has to be unpacked
   into a thousand files.
2. **Lossless structured fields.** :class:`~grkit.reader.Chapter` exposes metadata,
   sections, references, and cross-references, and :meth:`Chapter.to_record` emits a
   JSON-ready document.
3. **Faithful tables.** :class:`~grkit.model.Table` expands ``rowspan``/``colspan`` into a
   real grid, keeps multi-row headers, and renders back to HTML with the original layout.
   :mod:`grkit.testing` builds the "testing strategy vs yield" view on top of it --
   sequencing versus copy-number testing -- while preserving the source table for display.

Quick start::

    from grkit import Corpus

    with Corpus("gene_NBK1116.tar.gz") as corpus:      # or a directory of .nxml
        print(len(corpus))                             # 996
        chapter = corpus.get("abs")
        print(chapter.genes, chapter.title)
        for result in extract_testing_yield(chapter):
            print(result.comparison())                 # sequencing vs CNV per gene
            print(result.to_html())                    # faithful HTML table
"""

from __future__ import annotations

from .corpus import Corpus
from .inline import (
    InlineOptions,
    element_to_html,
    element_to_text,
    element_to_text_and_marks,
)
from .model import (
    Author,
    Cell,
    CrossReference,
    Footnote,
    Link,
    Reference,
    Section,
    Table,
)
from .reader import Chapter, classify_url, parse_chapter, section_name
from .render import (
    DEFAULT_TABLE_CSS,
    footnotes_to_html,
    html_document,
    table_to_html,
    table_to_markdown,
)
from .sources import (
    ChapterSource,
    DirectorySource,
    FileSource,
    TarGzSource,
    is_archive,
    open_source,
)
from .tables import parse_table_wrap, table_to_record, tables_of
from .testing import (
    MethodKind,
    TestingYieldTable,
    YieldObservation,
    YieldRow,
    classify_method,
    csv_text,
    extract_all,
    extract_from_table,
    extract_testing_yield,
    is_testing_yield_table,
    parse_yield,
    to_records,
    write_csv,
    write_json,
)
from .updater import (
    ArchiveInfo,
    ChangeReport,
    fetch_archive,
    is_stale,
    resolve_archive,
    update_corpus,
)

__version__ = "0.1.0"

#: Convenience alias for the most common entry point.
chapter_testing = extract_testing_yield

__all__ = [
    "__version__",
    # corpus & sources
    "Corpus",
    "ChapterSource",
    "DirectorySource",
    "FileSource",
    "TarGzSource",
    "open_source",
    "is_archive",
    # chapter model
    "Chapter",
    "parse_chapter",
    "section_name",
    "classify_url",
    "Author",
    "CrossReference",
    "Footnote",
    "Link",
    "Reference",
    "Section",
    # tables
    "Table",
    "Cell",
    "parse_table_wrap",
    "tables_of",
    "table_to_record",
    # rendering
    "table_to_html",
    "table_to_markdown",
    "footnotes_to_html",
    "html_document",
    "DEFAULT_TABLE_CSS",
    "element_to_text",
    "element_to_html",
    "element_to_text_and_marks",
    "InlineOptions",
    # testing strategy vs yield
    "TestingYieldTable",
    "YieldRow",
    "YieldObservation",
    "MethodKind",
    "extract_testing_yield",
    "chapter_testing",
    "extract_from_table",
    "extract_all",
    "is_testing_yield_table",
    "classify_method",
    "parse_yield",
    "csv_text",
    "to_records",
    "write_json",
    "write_csv",
    # self-update
    "ArchiveInfo",
    "ChangeReport",
    "resolve_archive",
    "fetch_archive",
    "is_stale",
    "update_corpus",
]
