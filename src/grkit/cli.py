"""Command-line interface: ``grkit``.

Subcommands::

    grkit update                    # self-update: fetch/cache the NCBI archive
    grkit update --check            # report whether a newer archive exists
    grkit list                      # chapters with genes/diseases
    grkit show abs                  # full structured record for one chapter
    grkit tables -c abs             # every table, as JSON (or --html)
    grkit testing -c abs            # testing strategy vs yield
    grkit testing --strategy        # add the normalized sequencing-vs-CNV view
    grkit testing -f csv -o out.csv # flat records for spreadsheets
    grkit html -c abs -o abs.html   # standalone HTML with faithful tables

The corpus is located automatically: ``--source`` wins, then ``$GRKIT_SOURCE``, then a
cached ``NBK1116.tar.gz`` in ``--cache-dir``/``$GRKIT_CACHE`` (default
``~/.cache/grkit``), then a local ``nxmls/`` directory.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections.abc import Iterable
from typing import Any

from .corpus import Corpus
from .reader import Chapter
from .render import html_document
from .sources import open_source
from .updater import GENEREVIEWS_ACCESSION, fetch_archive, is_stale, resolve_archive

DEFAULT_CACHE_DIR = os.environ.get("GRKIT_CACHE") or os.path.join(
    os.path.expanduser("~"), ".cache", "grkit"
)


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────
def _emit(text: str, path: str | None) -> None:
    """Write ``text`` to ``path``, or to stdout when no path is given."""
    if path:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        print(f"wrote {path}", file=sys.stderr)
    else:
        sys.stdout.write(text)
        if not text.endswith("\n"):
            sys.stdout.write("\n")


def _dumps(value: Any) -> str:
    return json.dumps(value, indent=2, ensure_ascii=False)


def _cached_archive(cache_dir: str) -> str:
    return os.path.join(cache_dir, f"{GENEREVIEWS_ACCESSION}.tar.gz")


def resolve_source(args: argparse.Namespace, *, allow_download: bool = False) -> str:
    """Pick the corpus location from CLI flags, env vars, or the cache."""
    if getattr(args, "source", None):
        return args.source
    env = os.environ.get("GRKIT_SOURCE")
    if env:
        return env
    cached = _cached_archive(getattr(args, "cache_dir", DEFAULT_CACHE_DIR))
    if os.path.exists(cached):
        return cached
    if allow_download:
        info = fetch_archive(getattr(args, "cache_dir", DEFAULT_CACHE_DIR))
        return info.local_path or cached
    for candidate in ("nxmls", "gene_NBK1116"):
        if os.path.isdir(candidate):
            return candidate
    raise SystemExit(
        "No corpus found. Run `grkit update` first, or pass --source "
        "(a .tar.gz archive, an .nxml directory, or a single .nxml file)."
    )


def _chapters(corpus: Corpus, args: argparse.Namespace) -> Iterable[Chapter]:
    """Iterate chapters honouring ``--chapter`` and ``--limit``."""
    if getattr(args, "chapter", None):
        for chapter_id in args.chapter:
            chapter = corpus.get(chapter_id)
            if chapter is None:
                print(f"warning: chapter {chapter_id!r} not found", file=sys.stderr)
                continue
            yield chapter
        return
    limit = getattr(args, "limit", None)
    count = 0
    for chapter in corpus:
        yield chapter
        count += 1
        if limit is not None and count >= limit:
            return


# ──────────────────────────────────────────────────────────────────────────────
# Commands
# ──────────────────────────────────────────────────────────────────────────────
def cmd_update(args: argparse.Namespace) -> int:
    if args.check:
        stale = is_stale(args.cache_dir, accession=args.accession)
        info = resolve_archive(args.accession)
        print(
            json.dumps(
                {
                    "accession": args.accession,
                    "remote_last_updated": info.last_updated,
                    "remote_url": info.url,
                    "cache_dir": args.cache_dir,
                    "cached": os.path.exists(_cached_archive(args.cache_dir)),
                    "stale": stale,
                },
                indent=2,
            )
        )
        return 0

    info = fetch_archive(
        args.cache_dir, accession=args.accession, force=args.force, quiet=args.quiet
    )
    print(json.dumps(info.__dict__, indent=2))
    return 0


def cmd_list(args: argparse.Namespace) -> int:
    corpus = Corpus(resolve_source(args))
    records = []
    for chapter in _chapters(corpus, args):
        records.append(
            {
                "id": chapter.id,
                "title": chapter.title,
                "genes": chapter.genes,
                "diseases": chapter.diseases,
                "is_retired": chapter.is_retired,
                "content_type": chapter.content_type,
            }
        )
    _emit(_dumps(records) if args.json else "\n".join(
        f"{r['id']}\t{r['title']}" for r in records
    ), args.output)
    return 0


def cmd_show(args: argparse.Namespace) -> int:
    corpus = Corpus(resolve_source(args))
    records = [chapter.to_record(table_html=args.html) for chapter in _chapters(corpus, args)]
    payload: Any = records[0] if len(records) == 1 else records
    _emit(_dumps(payload), args.output)
    return 0


def cmd_tables(args: argparse.Namespace) -> int:
    from .tables import table_to_record

    source = resolve_source(args)
    with open_source(source) as src:
        corpus = Corpus(src)
        records = []
        for chapter in _chapters(corpus, args):
            for table in chapter.tables:
                if args.testing_only:
                    from .testing import is_testing_yield_table

                    if not is_testing_yield_table(table):
                        continue
                record = table_to_record(table, include_html=args.html)
                record["chapter_title"] = chapter.title
                records.append(record)

    if args.format == "html":
        body = "\n".join(r.get("html", "") for r in records)
        _emit(html_document(body, title="GeneReviews tables"), args.output)
    else:
        _emit(_dumps(records), args.output)
    return 0


def cmd_testing(args: argparse.Namespace) -> int:
    from .testing import csv_text, extract_testing_yield, to_records, write_csv

    source = resolve_source(args)
    results = []
    with open_source(source) as src:
        corpus = Corpus(src)
        for chapter in _chapters(corpus, args):
            for result in extract_testing_yield(
                chapter,
                require_caption_match=args.require_caption,
                min_confidence=args.min_confidence,
            ):
                results.append(result)

    if args.format == "csv":
        if args.output:
            write_csv(results, args.output)
            print(f"wrote {args.output}", file=sys.stderr)
        else:
            sys.stdout.write(csv_text(to_records(results)))
        return 0

    if args.format == "html":
        parts: list[str] = []
        for result in results:
            parts.append(
                f"<h2>{result.chapter_id}: {result.caption or result.table_id}</h2>"
            )
            if args.strategy:
                parts.append("<h3>Strategy vs yield</h3>")
                parts.append(result.comparison_to_html())
            parts.append("<h3>Source table</h3>")
            parts.append(result.to_html())
        page = html_document(
            "\n".join(parts), title="GeneReviews testing strategy vs yield"
        )
        _emit(page, args.output)
        return 0

    payload = [r.to_dict(include_source_html=args.html) for r in results]
    _emit(_dumps(payload), args.output)
    return 0


def cmd_html(args: argparse.Namespace) -> int:
    from .testing import extract_testing_yield

    source = resolve_source(args)
    parts: list[str] = []
    with open_source(source) as src:
        corpus = Corpus(src)
        for chapter in _chapters(corpus, args):
            parts.append(f'<h1 id="{chapter.id}">{chapter.title}</h1>')
            if chapter.genes:
                parts.append(f'<p class="gr-genes">Genes: {", ".join(chapter.genes)}</p>')
            results = extract_testing_yield(chapter)
            if not results:
                parts.append("<p><em>No testing-strategy table found.</em></p>")
                continue
            for result in results:
                parts.append(f"<h2>{result.label} {result.caption}</h2>")
                if args.strategy:
                    parts.append(result.comparison_to_html())
                parts.append(result.to_html())
    _emit(html_document("\n".join(parts), title="GeneReviews"), args.output)
    return 0


# ──────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ──────────────────────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="grkit",
        description=(
            "Read the GeneReviews NXML corpus "
            "(structured fields + testing-strategy tables)."
        ),
    )
    parser.add_argument("--version", action="store_true", help="print version and exit")

    # Source selection is accepted after any subcommand.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--source",
        help="corpus location: a .tar.gz archive, an .nxml directory, or a single .nxml file",
    )
    common.add_argument(
        "--cache-dir", default=DEFAULT_CACHE_DIR, help="cache dir for the NCBI archive"
    )

    sub = parser.add_subparsers(dest="command", required=False)

    def add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("-c", "--chapter", action="append", help="chapter id (repeatable)")
        sp.add_argument("--limit", type=int, help="stop after N chapters")
        sp.add_argument("-o", "--output", help="write output to this file")

    p_update = sub.add_parser("update", help="fetch/cache the archive from NCBI", parents=[common])
    p_update.add_argument("--accession", default=GENEREVIEWS_ACCESSION)
    p_update.add_argument("--force", action="store_true", help="re-download unconditionally")
    p_update.add_argument("--check", action="store_true", help="only report staleness")
    p_update.add_argument("-q", "--quiet", action="store_true")
    p_update.set_defaults(func=cmd_update)

    p_list = sub.add_parser("list", help="list chapters", parents=[common])
    add_common(p_list)
    p_list.add_argument("--json", action="store_true", help="JSON instead of TSV")
    p_list.set_defaults(func=cmd_list)

    p_show = sub.add_parser("show", help="print a chapter's structured record", parents=[common])
    add_common(p_show)
    p_show.add_argument("--html", action="store_true", help="embed table HTML")
    p_show.set_defaults(func=cmd_show)

    p_tables = sub.add_parser("tables", help="print tables", parents=[common])
    add_common(p_tables)
    p_tables.add_argument("--testing-only", action="store_true", help="only testing tables")
    p_tables.add_argument("--html", action="store_true", help="embed faithful table HTML")
    p_tables.add_argument("-f", "--format", choices=("json", "html"), default="json")
    p_tables.set_defaults(func=cmd_tables)

    p_testing = sub.add_parser(
        "testing", help="testing strategy vs yield", parents=[common]
    )
    add_common(p_testing)
    p_testing.add_argument("-f", "--format", choices=("json", "csv", "html"), default="json")
    p_testing.add_argument(
        "--strategy", action="store_true", help="include the sequencing-vs-CNV view"
    )
    p_testing.add_argument("--html", action="store_true", help="embed source table HTML (json)")
    p_testing.add_argument(
        "--require-caption", action="store_true", help="only tables whose caption names testing"
    )
    p_testing.add_argument("--min-confidence", type=float, default=0.0)
    p_testing.set_defaults(func=cmd_testing)

    p_html = sub.add_parser(
        "html", help="render chapters' testing tables to one HTML page", parents=[common]
    )
    add_common(p_html)
    p_html.add_argument("--strategy", action="store_true", help="also show normalized view")
    p_html.set_defaults(func=cmd_html)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if getattr(args, "version", False):
        from . import __version__

        print(f"grkit {__version__}")
        return 0

    if not getattr(args, "command", None):
        parser.print_help()
        return 0

    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
