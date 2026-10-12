"""Shared pytest fixtures for the grkit test suite.

The tests are fully hermetic: they use the committed synthetic fixture
``fixtures/sample-chapter.nxml`` and never touch the network or the real corpus.
"""

from __future__ import annotations

import pathlib
import sys

import pytest

# Allow running the suite from a checkout without installing the package.
SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

FIXTURES = pathlib.Path(__file__).resolve().parent / "fixtures"
SAMPLE = FIXTURES / "sample-chapter.nxml"


@pytest.fixture(scope="session")
def fixture_dir() -> pathlib.Path:
    return FIXTURES


@pytest.fixture(scope="session")
def sample_path() -> pathlib.Path:
    return SAMPLE


@pytest.fixture
def chapter():
    """The parsed fixture chapter (``id='tst'``)."""
    from grkit.reader import Chapter

    assert SAMPLE.exists(), f"missing fixture: {SAMPLE}"
    return Chapter.from_file(str(SAMPLE))


@pytest.fixture
def table_of(chapter):
    """Look up a table by id, failing loudly when absent."""

    def _get(table_id: str):
        for table in chapter.tables:
            if table.id == table_id:
                return table
        raise AssertionError(f"table {table_id!r} not found; have {[t.id for t in chapter.tables]}")

    return _get
