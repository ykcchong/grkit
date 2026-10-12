"""Read a GeneReviews NXML chapter into structured, JSON-serializable fields.

A chapter is a ``<book-part-wrapper>`` (JATS/Book DTD v2.0) whose metadata lives in
``<book-part-meta>`` and whose content lives in ``<body>``/``<back>``. This module maps
that onto a small, stable API:

    >>> ch = Chapter.from_file("nxmls/abs.nxml")
    >>> ch.title
    'Cytochrome P450 Oxidoreductase Deficiency'
    >>> ch.genes
    ['POR']
    >>> ch.section_names
    ['GeneReview_Scope', 'Diagnosis', ...]
    >>> diagnosis = ch.get_section("Diagnosis")
    >>> diagnosis.tables[0].to_dicts()

Keyword groups carry the indexed vocabulary (``GeneSymbol``, ``DiseaseName``, ``Include``,
``Synonym``, ``ProteinName``), and are read case-insensitively so small DTD variations
between chapters do not silently drop metadata.
"""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from typing import Any

from .inline import element_to_html, element_to_text, localname
from .model import Author, CrossReference, Link, Reference, Section, Table
from .tables import tables_of

_XLINK_HREF = "{http://www.w3.org/1999/xlink}href"

# External-resource URL patterns, checked in order. The first group is the bare identifier.
_RESOURCE_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("omim", re.compile(r"omim\.org/entry/(\d+)", re.I)),
    ("clinvar", re.compile(r"clinvar/(?:variation/)?(\d+)", re.I)),
    ("gtr", re.compile(r"gtr/(?:tests/)?(\d+)", re.I)),
    ("genereviews", re.compile(r"genereviews/(?:NBK\d+/?)?([A-Za-z0-9\-_]+)", re.I)),
    ("pmid", re.compile(r"pubmed(?:\.ncbi\.nlm\.nih\.gov)?/(\d+)", re.I)),
    ("hgnc", re.compile(r"hgnc[^/]*/(?:gene/)?([A-Za-z0-9\-_]+)", re.I)),
    ("medgen", re.compile(r"medgen/(\d+)", re.I)),
    (
        "orphanet",
        re.compile(r"orpha\.net/(?:consor/cgi-bin/OC_Exp\.php\?lng=EN&Expert=)?(\d+)", re.I),
    ),
    ("gene", re.compile(r"ncbi\.nlm\.nih\.gov/gene/(\d+)", re.I)),
    ("doi", re.compile(r"doi\.org/(10\.\S+)", re.I)),
)


def section_name(section_id: str) -> str:
    """``"abs.Diagnosis"`` -> ``"Diagnosis"``."""
    return section_id.rsplit(".", 1)[-1] if "." in section_id else section_id


def sections_of(
    scope: ET.Element | None, *, depth: int = 0, chapter_id: str = ""
) -> list[Section]:
    """Return the direct ``<sec>`` children of ``scope`` as :class:`Section` objects."""
    if scope is None:
        return []
    out: list[Section] = []
    for sec in scope.findall("sec"):
        sec_id = sec.get("id", "")
        out.append(
            Section(
                id=sec_id,
                name=section_name(sec_id),
                title=element_to_text(sec.find("title")),
                depth=depth,
                element=sec,
                chapter_id=chapter_id,
            )
        )
    return out


def classify_url(url: str) -> tuple[str, str]:
    """Return ``(db, value)`` for an external URL, or ``("uri", url)`` if unknown."""
    if not url:
        return "", ""
    for db, pattern in _RESOURCE_PATTERNS:
        match = pattern.search(url)
        if match:
            return db, match.group(1)
    return "uri", url


class Chapter:
    """A single GeneReviews chapter, parsed from its NXML."""

    def __init__(self, root: ET.Element, *, source_path: str | None = None) -> None:
        self._root = root
        self.source_path = source_path

    # ── construction ─────────────────────────────────────────────────────────
    @classmethod
    def from_file(cls, filepath: str) -> Chapter:
        """Parse an ``.nxml`` file (or an open binary file object)."""
        if hasattr(filepath, "read"):
            root = ET.parse(filepath).getroot()  # type: ignore[arg-type]
            return cls(root)
        root = ET.parse(filepath).getroot()
        return cls(root, source_path=str(filepath))

    @classmethod
    def from_string(cls, xml_text: str | bytes) -> Chapter:
        return cls(ET.fromstring(xml_text))

    @classmethod
    def try_from_file(cls, filepath: str) -> Chapter | None:
        """Parse ``filepath``, returning ``None`` instead of raising on malformed XML."""
        try:
            return cls.from_file(filepath)
        except (ET.ParseError, OSError, ValueError):
            return None

    # ── low-level access ─────────────────────────────────────────────────────
    @property
    def root(self) -> ET.Element:
        return self._root

    @property
    def meta(self) -> ET.Element | None:
        return self._root.find(".//book-part-meta")

    @property
    def body(self) -> ET.Element | None:
        return self._root.find(".//body")

    @property
    def back(self) -> ET.Element | None:
        return self._root.find(".//back")

    def _keyword_group(self, *group_types: str) -> list[str]:
        """Values of ``<kwd>`` in every matching ``kwd-group`` (case-insensitive)."""
        wanted = {t.lower() for t in group_types}
        values: list[str] = []
        for group in self._root.iter("kwd-group"):
            group_type = (group.get("kwd-group-type") or "").lower()
            if group_type not in wanted:
                continue
            for kwd in group.findall("kwd"):
                text = element_to_text(kwd)
                if text and text not in values:
                    values.append(text)
        return values

    def _related_objects(self, link_type: str) -> list[dict[str, str]]:
        out: list[dict[str, str]] = []
        for obj in self._root.iter("related-object"):
            if (obj.get("link-type") or "") != link_type:
                continue
            out.append(
                {
                    "document_id": obj.get("document-id", ""),
                    "document_type": obj.get("document-type", ""),
                    "title": element_to_text(obj),
                }
            )
        return out

    # ── identity & metadata ──────────────────────────────────────────────────
    @property
    def id(self) -> str:
        """Chapter id, e.g. ``"abs"`` (from the wrapper id, else the filename stem)."""
        return self._root.get("id") or (
            self.source_path.rsplit("/", 1)[-1].removesuffix(".nxml") if self.source_path else ""
        )

    @property
    def content_type(self) -> str:
        return self._root.get("content-type", "")

    @property
    def title(self) -> str:
        title = self._root.find(".//title-group/title")
        if title is not None:
            return element_to_text(title)
        return element_to_text(self._root.find(".//book-part-meta//title"))

    @property
    def alt_titles(self) -> list[str]:
        return [
            element_to_text(el)
            for el in self._root.findall(".//title-group/alt-title")
            if element_to_text(el)
        ]

    @property
    def synonyms(self) -> list[str]:
        return self._keyword_group("synonym", "synonyms")

    @property
    def genes(self) -> list[str]:
        return self._keyword_group("genesymbol", "gene")

    @property
    def diseases(self) -> list[str]:
        """Disease names the chapter is indexed under."""
        return self._keyword_group("diseasename", "disease")

    @property
    def included_phenotypes(self) -> list[str]:
        return self._keyword_group("include")

    @property
    def proteins(self) -> list[str]:
        return self._keyword_group("proteinname")

    @property
    def authors(self) -> list[Author]:
        out: list[Author] = []
        meta = self.meta
        if meta is None:
            return out
        for contrib in meta.findall(".//contrib-group/contrib"):
            name = contrib.find("name")
            if name is None:
                continue
            email = contrib.findtext("email") or ""
            out.append(
                Author(
                    surname=element_to_text(name.find("surname")),
                    given_names=element_to_text(name.find("given-names")),
                    email=email.strip(),
                    affiliation=element_to_text(contrib.find("aff")),
                )
            )
        return out

    def _date(self, date_type: str) -> str:
        for date in self._root.iter("date"):
            if (date.get("date-type") or "") != date_type:
                continue
            year = date.findtext("year") or ""
            month = (date.findtext("month") or "").zfill(2)
            day = (date.findtext("day") or "").zfill(2)
            parts = [p for p in (year, month, day) if p and p.strip("0")]
            return "-".join(parts)
        return ""

    @property
    def created_date(self) -> str:
        return self._date("created")

    @property
    def updated_date(self) -> str:
        return self._date("updated")

    @property
    def is_retired(self) -> bool:
        """GeneReviews marks retirement in the chapter's own title."""
        return "RETIRED CHAPTER" in self.title.upper()

    @property
    def previous_chapter(self) -> list[dict[str, str]]:
        return self._related_objects("previous-part-link")

    @property
    def next_chapter(self) -> list[dict[str, str]]:
        return self._related_objects("next-part-link")

    # ── content ──────────────────────────────────────────────────────────────
    @property
    def sections(self) -> list[Section]:
        """Top-level body sections in document order."""
        return sections_of(self.body, depth=0, chapter_id=self.id)

    @property
    def section_names(self) -> list[str]:
        return [s.name for s in self.sections]

    def get_section(self, name: str) -> Section | None:
        """Find a section by short name (``"Diagnosis"``) or full id, case-insensitively."""
        wanted = name.lower()
        for sec in self.sections:
            if sec.name.lower() == wanted or sec.id.lower() == wanted:
                return sec
        # Fall back to a depth-first search of nested sections.
        for sec in self.sections:
            for child in self._walk(sec):
                if child.name.lower() == wanted or child.id.lower() == wanted:
                    return child
        return None

    def _walk(self, section: Section) -> Iterator[Section]:
        for child in section.subsections:
            yield child
            yield from self._walk(child)

    @property
    def summary(self) -> list[dict[str, str]]:
        """The structured abstract: one entry per abstract subsection."""
        abstract = self._root.find(".//book-part-meta/abstract")
        if abstract is None:
            abstract = self._root.find(".//abstract")
        if abstract is None:
            return []
        out: list[dict[str, str]] = []
        default_title = element_to_text(abstract.find("title"))
        for sec in abstract.findall("sec"):
            out.append(
                {
                    "title": element_to_text(sec.find("title")),
                    "content": element_to_text(sec),
                    "html": element_to_html(sec),
                }
            )
        if not out:
            content = element_to_text(abstract)
            if default_title and content.startswith(default_title):
                content = content[len(default_title) :].strip()
            if content:
                out.append(
                    {
                        "title": default_title,
                        "content": content,
                        "html": element_to_html(abstract),
                    }
                )
        return out

    @property
    def references(self) -> list[Reference]:
        out: list[Reference] = []
        for ref in self._root.iter("ref"):
            citation = ref.find("element-citation")
            if citation is None:
                citation = ref.find("mixed-citation")
            if citation is None:
                continue
            authors: list[str] = []
            for name in citation.iter("name"):
                surname = element_to_text(name.find("surname"))
                given = element_to_text(name.find("given-names"))
                if surname:
                    authors.append(f"{surname} {given}".strip())
            pmid = ""
            for pub_id in citation.iter("pub-id"):
                if (pub_id.get("pub-id-type") or "").lower() == "pmid":
                    pmid = element_to_text(pub_id)
                    break
            out.append(
                Reference(
                    id=ref.get("id", ""),
                    authors=authors,
                    title=element_to_text(citation.find("article-title")),
                    source=element_to_text(citation.find("source")),
                    year=element_to_text(citation.find("year")),
                    volume=element_to_text(citation.find("volume")),
                    fpage=element_to_text(citation.find("fpage")),
                    lpage=element_to_text(citation.find("lpage")),
                    pmid=pmid,
                    doi=next(
                        (
                            element_to_text(p)
                            for p in citation.iter("pub-id")
                            if (p.get("pub-id-type") or "").lower() == "doi"
                        ),
                        "",
                    ),
                    pub_type=citation.get("publication-type", ""),
                )
            )
        return out

    @property
    def links(self) -> list[Link]:
        """Every ``<ext-link>``/``<uri>`` in the document, with its resolved href."""
        out: list[Link] = []
        for el in self._root.iter():
            name = localname(el.tag)
            if name not in ("ext-link", "uri"):
                continue
            href = el.get(_XLINK_HREF) or el.get("href") or ""
            if not href:
                continue
            out.append(
                Link(
                    text=element_to_text(el),
                    href=href,
                    link_type=el.get("ext-link-type", "") or el.get("link-type", ""),
                )
            )
        return out

    @property
    def cross_references(self) -> list[CrossReference]:
        """Deduplicated links classified into external databases (OMIM, ClinVar, GTR, ...)."""
        seen: set[tuple[str, str]] = set()
        out: list[CrossReference] = []
        for link in self.links:
            db, value = classify_url(link.href)
            if not db or (db, value) in seen:
                continue
            seen.add((db, value))
            out.append(CrossReference(db=db, value=value, text=link.text, url=link.href))
        return out

    # ── tables ───────────────────────────────────────────────────────────────
    @property
    def tables(self) -> list[Table]:
        """Every table in the body, tagged with the top-level section that owns it."""
        out: list[Table] = []
        for sec in self.sections:
            out.extend(tables_of(sec.element, chapter_id=self.id, section=sec.name))
        return out

    def table_by_id(self, table_id: str) -> Table | None:
        for table in self.tables:
            if table.id == table_id:
                return table
        return None

    # ── structured output ────────────────────────────────────────────────────
    @property
    def body_text(self) -> str:
        """Full whitespace-collapsed plain text of the body, in document order."""
        return element_to_text(self.body) if self.body is not None else ""

    @property
    def text_length(self) -> int:
        return len(self.body_text)

    def section_tree(self) -> list[dict[str, Any]]:
        """Nested section dicts: id, name, title, depth, paragraphs, subsections."""

        def node(section: Section) -> dict[str, Any]:
            return {
                "id": section.id,
                "name": section.name,
                "title": section.title,
                "depth": section.depth,
                "paragraphs": section.paragraphs,
                "subsections": [node(child) for child in section.subsections],
            }

        return [node(sec) for sec in self.sections]

    def to_record(
        self, *, include_tables: bool = True, table_html: bool = True
    ) -> dict[str, Any]:
        """A complete, JSON-serializable representation of the chapter."""
        from .tables import table_to_record

        record: dict[str, Any] = {
            "id": self.id,
            "content_type": self.content_type,
            "file": (self.source_path or "").rsplit("/", 1)[-1],
            "title": self.title,
            "alt_titles": self.alt_titles,
            "synonyms": self.synonyms,
            "genes": self.genes,
            "diseases": self.diseases,
            "included_phenotypes": self.included_phenotypes,
            "proteins": self.proteins,
            "authors": [
                {
                    "surname": a.surname,
                    "given_names": a.given_names,
                    "full_name": a.full_name,
                    "email": a.email,
                }
                for a in self.authors
            ],
            "created_date": self.created_date,
            "updated_date": self.updated_date,
            "is_retired": self.is_retired,
            "text_length": self.text_length,
            "summary": self.summary,
            "section_names": self.section_names,
            "sections": self.section_tree(),
            "references": [
                {
                    "id": r.id,
                    "authors": r.authors,
                    "title": r.title,
                    "source": r.source,
                    "year": r.year,
                    "volume": r.volume,
                    "fpage": r.fpage,
                    "lpage": r.lpage,
                    "pmid": r.pmid,
                    "doi": r.doi,
                }
                for r in self.references
            ],
            "cross_references": [
                {"db": x.db, "value": x.value, "text": x.text, "url": x.url}
                for x in self.cross_references
            ],
        }
        if include_tables:
            record["tables"] = [
                table_to_record(t, include_html=table_html) for t in self.tables
            ]
        return record

    def to_json(self, *, indent: int | None = 2, **kwargs: Any) -> str:
        return json.dumps(self.to_record(**kwargs), indent=indent, ensure_ascii=False)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"<Chapter {self.id!r} {self.title!r} sections={len(self.sections)}>"


def parse_chapter(filepath: str) -> Chapter:
    """Convenience wrapper around :meth:`Chapter.from_file`."""
    return Chapter.from_file(filepath)
