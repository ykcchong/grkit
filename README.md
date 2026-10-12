# grkit

A dependency-free Python library for the [GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/)
NXML corpus (NCBI Bookshelf accession `NBK1116`).

It does four things:

1. **Updates itself.** Resolves the current archive URL from NCBI, caches
   `gene_NBK1116.tar.gz`, and re-downloads only when NCBI publishes a new version.
2. **Reads the archive in place.** Chapters are parsed *directly out of the `.tar.gz`* —
   no need to unpack ~1,000 files and ~870 MB onto disk — or from an extracted directory.
3. **Exposes structured fields.** Metadata, sections, tables, references, and
   cross-references, with a JSON-ready `to_record()`.
4. **Preserves tables losslessly, especially "testing strategy vs yield."** Real
   `rowspan`/`colspan` geometry, multi-row headers, and footnotes survive, and render back
   to faithful HTML — plus a normalized **sequencing vs copy-number** view.

Everything is stdlib-only (`xml.etree`, `tarfile`, `urllib`); `pytest` and `ruff` are the
only dev dependencies.

---

## Install

```bash
cd grkit
pip install -e .          # or: pip install -e ".[dev]"
```

Or run straight from a checkout with `PYTHONPATH=src python -m grkit.cli ...`.

## 1. Self-update

```bash
grkit update             # fetch/cache the archive (skips the download if unchanged)
grkit update --check     # report the remote version and whether the cache is stale
grkit update --force     # re-download unconditionally
```

```python
from grkit import fetch_archive, is_stale

info = fetch_archive("~/.cache/grkit")      # ArchiveInfo with .local_path, .last_updated
print(info.url, info.last_updated)
if is_stale("~/.cache/grkit"):              # cheap, offline-safe staleness check
    fetch_archive("~/.cache/grkit", force=True)
```

The FTP path is content-hashed (`.../ca/84/gene_NBK1116.tar.gz`) and changes whenever NCBI
re-shards the archive, so `resolve_archive()` reads the published `file_list.csv` index
keyed by the stable accession `NBK1116` instead of hard-coding a URL. `is_stale()` fails
open (returns `True`) when offline, so callers refresh rather than silently serving stale
data. A SHA-256 manifest and add/modify/remove diff are also available via `update_corpus`
for pipelines that do want an extracted working tree.

## 2. Read chapters from the tarball or a directory

```python
from grkit import Corpus, open_source

with Corpus("gene_NBK1116.tar.gz") as corpus:     # ~620 MB tarball, read in place
    len(corpus)                                   # 996
    chapter = corpus.get("abs")
    chapter.title, chapter.genes
    for ch in corpus.iter_chapters():             # single streaming pass, no temp files
        ...
    corpus.find_by_gene("POLG")
    corpus.search("riboflavin")

with open_source("nxmls/") as source:             # or an extracted directory
    source.names()[:3]
```

Two access patterns are offered because gzip is not seekable:

| Pattern | Cost | Use for |
|---|---|---|
| `iter_chapters()` / `iter_bytes()` | one decompression pass, no extra disk | whole-corpus scans |
| `get("abs")` / `read("abs.nxml")` | a few rescans, then a temp spool | random lookup by id |

`get()` serves the first few lookups by rescanning (so a single `get("abs")` costs ~3 s and
writes nothing); after that it transparently decompresses once to a temporary tar for true
seeks, deleted on `close()`. `TarGzSource(..., spool=False)` never writes to disk.

## 3. Structured fields

```python
ch = corpus.get("abs")
ch.id, ch.title, ch.genes, ch.diseases, ch.synonyms, ch.proteins
ch.authors, ch.created_date, ch.updated_date, ch.is_retired
ch.section_names                       # ['GeneReview_Scope', 'Diagnosis', ...]
ch.get_section("Diagnosis").paragraphs
ch.references[0].pmid
[(x.db, x.value) for x in ch.cross_references]   # [('omim', '124015'), ('clinvar', '...')]
ch.summary                             # structured abstract

record = ch.to_record()                # JSON-serializable, includes tables as dicts
```

## 4. Tables: preserved, not flattened

The obvious way to read a table — walk `<tr>` and `<td>` in order — is **wrong** for this
corpus and quietly corrupts the exact table that matters most. In 697 of the 868 canonical
testing tables a gene cell uses `rowspan`, so a positional reader returns a short second
row and shifts the yield into the `Method` column:

```
source:   Gene | Method                                | Proportion ... Detectable by Method
          POR  | Sequence analysis                     | 92%
               | Gene-targeted deletion/duplication... | 2.5%     <- POR spans both rows

naive:    ['Gene-targeted deletion/duplication analysis', '2.5%']      # 2.5% is now "Method"
grkit:    ['POR', 'Gene-targeted deletion/duplication analysis', '2.5%']
```

`grkit` expands spans into a real rectangular grid, so alignment is never lost:

```python
table = ch.table_by_id("abs.T.molecular_genetic_testing_used_in")
table.to_grid(filled=True)     # rectangular, span values filled
table.to_dicts()               # keyed by header path
table.column_header_path(2)    # ['Proportion of ... Detectable by Method']
table.to_html()                # faithful HTML: real rowspan/colspan, scopes, footnotes
table.to_markdown()            # composite keys, since Markdown cannot merge cells
```

Three further details are preserved that flat readers drop:

- **Multi-row headers.** 172 tables put *methods in columns*
  (`... Identified by Method` over `Sequence analysis` / `Gene-targeted deletion/duplication
  analysis`, sometimes nested again as `Common deletions` / `Other deletions`).
  `column_header_path()` returns the full hierarchy.
- **Footnote markers.** `100% 4, 5` becomes text `100%` with `cell.marks == ('4', '5')`, and
  the marker renders as a link to the footnote. A `<sup>` is only treated as a marker when
  its value is a *defined* footnote label **and** it is detached from the preceding word —
  so `m<sup>2</sup>` stays a superscript and a cell that is just `1` is never emptied.
- **Inline markup.** `<italic>` gene names, `<break/>` line breaks, and links survive into
  the HTML instead of collapsing to flat text (which otherwise produces run-together
  strings like `...syndromeCongenital adrenal hyperplasia...`).

## 5. Testing strategy vs yield (sequencing vs copy number)

```python
from grkit import extract_testing_yield

for result in extract_testing_yield(ch):
    result.layout            # 'long' (Method column) or 'wide' (methods are columns)
    result.confidence        # how sure we are this is the right table
    result.methods           # ['Sequence analysis', 'Gene-targeted deletion/duplication analysis']
    for row in result.rows:
        row.gene, row.method, row.method_kind, row.yield_raw
        row.is_sequencing, row.is_copy_number
        row.observation.percent, row.observation.low, row.observation.high

    result.comparison()      # per gene: sequencing yield vs copy-number yield
    result.to_html()         # the ORIGINAL table, layout-faithful
    result.comparison_to_html()
```

```python
>>> result.comparison()
[{'gene': 'POR',
  'sequencing': '92%',   'sequencing_percent': 92.0,
  'copy_number': '2.5%', 'copy_number_percent': 2.5,
  'other_methods': [], 'attribution': ''}, ...]
```

Methods are normalized to `sequence`, `deletion_duplication`, `cytogenetic`, `methylation`,
`repeat_expansion`, `panel_or_genome`, `targeted`, or `other`. Yields are parsed rather than
assumed to be percentages — `100%`, `<1%`, `>99%`, `~85%`, `5-15%`, `70%-80%`, `15/15`,
`Up to 100%`, `Nearly 100%`, `See footnote 7.`, `1 reported`, `NA`, and `Unknown` are all
distinguished, with `low`/`high` bounds and a `midpoint` for comparison.

Because the source `Table` is retained on every result, the published table can always be
rendered exactly as-is alongside the normalized view.

## 6. Gene → testing-coverage tables

Project the whole corpus onto a single gene: what testing is recommended, and what is the
expected yield?

```python
with Corpus("gene_NBK1116.tar.gz") as corpus:
    cov = corpus.gene_coverage("SCN1A")

cov.gene_ids                # ['6323']  -- from the entry's Table A gene link
cov.summary_rows()          # one row per method, with the yield range across entries
cov.sensitivity()           # {'sequencing': {...}, 'copy_number': {...}, 'preferred': 'sequence'}
cov.entries_with_table      # [fhm, gefs]
cov.entries_indexed_only    # entries that discuss the gene but publish no table
cov.first_line()            # the first-listed method, i.e. the preferred strategy
cov.coverage_matrix()       # {chapter_id: {method: yield}}
cov.to_html()               # coverage table (summary + per-entry detail)
```

```
SCN1A (NCBI Gene ID 6323)
  Sequence analysis                          sequence             entries=2  yields=['100%', '73%-92%']      82.5–100%
  Gene-targeted deletion/duplication analysis deletion_duplication entries=2  yields=['None reported', '8%-27%'] 0–17.5%
```

Four things make this more than a filter:

- **Gene cells are not single symbols.** They can be `HBA1 & HBA2`, a fifteen-symbol
  mitochondrial run, `Unknown`, or a whole sentence
  (`~0.8-Mb heterozygous deletion at 1q21.1`). `cell_genes()` splits and validates them, so
  querying `HBA2` finds the `HBA1 & HBA2` row. Source casing is used to reject prose, and
  the `C11orf65` family is handled as a legitimate exception.
- **Implicit attribution.** A single-gene entry such as `gefs` (SCN1A Seizure Disorders)
  has a testing table with no usable Gene column; those rows are attributed to the gene and
  flagged `implicit=True` rather than dropped (or wrongly spread over an overview's genes).
- **Coverage gaps are reported.** SCN5A is indexed in three entries but only `brugada`
  publishes a method/yield table — `dcm-ov` and `rws` come back from
  `entries_indexed_only`. "Discusses the gene, publishes no coverage data" is part of the
  answer.
- **Qualitative yields stay out of the numbers.** `None reported` is an explicit zero and
  joins the range; `1 family`, `Unknown` and `See footnote 7.` are surfaced separately so a
  range is never read as the whole story.

Query a panel in one pass with `corpus.gene_coverages(["SCN1A", "SCN5A", ...])` — scanning
dominates the cost (~6 s for 996 chapters), so 60 genes cost about the same as one.

```bash
grkit coverage SCN1A SCN5A HBA2 KCNQ1                     # JSON
grkit coverage SCN1A -f html --source-tables -o out.html  # coverage + original tables
grkit coverage SCN1A SCN5A -f csv -o coverage.csv
```

![Gene testing-coverage tables](out/coverage-preview.png)

## CLI

```bash
grkit update   [--check|--force]                 # self-update
grkit list     [--limit N] [--json]              # chapter ids and titles
grkit show     -c abs [--html]                   # full structured record
grkit tables   -c abs [--testing-only] [--html]  # tables (JSON or standalone HTML)
grkit testing  -c abs --strategy                 # testing strategy vs yield (JSON)
grkit testing  --strategy -f html -o out.html    # faithful tables + normalized view
grkit testing  -f csv -o tables.csv              # flat records
grkit html     -c abs -c a-thal --strategy -o page.html
grkit coverage SCN1A SCN5A [-f json|csv|html] [--source-tables] [-o out]
```

`--source` accepts a `.tar.gz`, an `.nxml` directory, or a single file; otherwise
`$GRKIT_SOURCE`, then the cache, then a local `nxmls/` directory is used.

## Generated examples (`out/`)

Produced from the real corpus with the commands above:

| File | What it is |
|---|---|
| [`out/testing-tables.html`](out/testing-tables.html) | Faithful HTML for four chapters (long, wide/nested, stray-header-rowspan, two-yield-column) plus the normalized view |
| [`out/testing-strategy-vs-yield.csv`](out/testing-strategy-vs-yield.csv) | 3,271 flat (gene, method, yield) observations across 864 tables |
| [`out/testing-strategy-vs-yield.json`](out/testing-strategy-vs-yield.json) | Same, with per-gene `strategy_vs_yield` comparisons |
| [`out/chapter-abs.json`](out/chapter-abs.json) | One chapter's full structured record, including table HTML |
| [`out/gene-testing-coverage.html`](out/gene-testing-coverage.html) | Gene → testing-coverage tables for SCN1A, SCN5A, HBA2, KCNQ1 |
| [`out/gene-coverage.csv`](out/gene-coverage.csv) | The same coverage rows in flat form |
| [`out/rendered-preview.png`](out/rendered-preview.png) | Screenshot of the rendered table HTML |
| [`out/coverage-preview.png`](out/coverage-preview.png) | Screenshot of the rendered coverage tables |

## Validation against the real corpus

```
996 chapters parsed                0 errors      5.2 s (from a directory)
864 testing tables extracted       863 chapters
716 long layout / 148 wide layout
859 of 864 tables confidence >= 0.8 (508 at 1.0)
method kinds: sequence 1470, deletion_duplication 1458, cytogenetic 160, other 79,
              targeted 68, repeat_expansion 17, methylation 14, panel_or_genome 5
yield kinds:  percent 1836, fraction 206, unknown 283, none_reported 360,
              not_applicable 186, count 105, see_footnote 91, text 154, multi_percent 13
```

Gene→coverage projection over a 60-gene panel: 60 coverage tables built in one 6.3 s scan,
0 exceptions; tokenization rejects prose, chromosome arms, karyotype complements and assay
abbreviations, leaving 1,267 distinct real gene symbols.

Chapters without a testing table are correctly excluded: appendices and front matter, plus
gene-attribution and differential-diagnosis tables that have a yield-looking column but no
method facet. Tables such as `Genes and Databases` (`Gene | Chromosome Locus | Protein`) are
deliberately rejected — `Chromosome Locus` must not be mistaken for a cytogenetic test.

## Tests and lint

```bash
python -m pytest tests/ -q     # 249 tests, hermetic (no network, no real corpus)
python -m ruff check src tests
```

Coverage includes span expansion and the row-misalignment regression, the stray single-row
`<thead>` `rowspan` quirk, footnote-marker detachment rules, HTML rendering (tag ordering,
escaping, anchor integrity), both testing-table layouts, yield/method parsing, tar streaming
and spooling, gene-cell tokenization and coverage projection, and the updater with a mocked network.

## Package layout

```
src/grkit/
  updater.py   NCBI URL resolution, cached download, manifest/diff, staleness
  sources.py   DirectorySource / TarGzSource / FileSource + open_source()
  corpus.py    Corpus: iterate, get, find_by_gene/disease, search
  reader.py    Chapter: metadata, sections, references, cross-references, to_record()
  model.py     Cell / Table / Footnote / Section / Reference / ...
  tables.py    span-aware <table-wrap> parsing, structured table records
  inline.py    rich inline -> plain text + HTML, footnote-marker handling
  render.py    faithful HTML / Markdown, footnote links, CSS
  testing.py   testing-strategy-vs-yield extraction and normalization
  coverage.py  gene -> testing-coverage tables (method/yield projection, gaps)
  cli.py       the grkit command
tests/         249 hermetic tests + synthetic fixtures covering every real structure
```

## Notes and limitations

- **Confidence is about table identity, not data quality.** An all-`Unknown` Table 1 scores
  high because it *is* the right table. Use `row.observation.kind` to judge yields.
- `rowspan` in a header is clamped so a header cell cannot spill into the body: GeneReviews
  emits stray `rowspan="2"` on single-row `<thead>` cells, and honouring it literally shifts
  every data column.
- Random access into a `.tar.gz` cannot be O(1) — the format is not seekable. Whole-corpus
  work should use `iter_chapters()`.
- Footnotes are linked with ids derived from the table id; pass `id_prefix=` when embedding
  several documents in one HTML page.
