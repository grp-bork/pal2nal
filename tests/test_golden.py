"""Compare the Python port against the golden outputs of pal2nal.v14.pl.

Every case in cases.tsv is run through the port and its stdout, stderr and
exit status are compared with what the Perl produced. The comparison is
byte for byte: the port is a reimplementation, not a rewrite, and any
difference in wording, spacing or stream is a regression.
"""

from __future__ import annotations

import io
import os
import shlex
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

TESTS = Path(__file__).parent
GOLDEN = TESTS / "golden"
EXPECTED = TESTS / "expected"


def load_cases() -> list[tuple[str, list[str]]]:
    cases = []
    for line in (TESTS / "cases.tsv").read_text().splitlines():
        if not line.strip():
            continue
        name, _, args = line.partition("\t")
        cases.append((name, shlex.split(args)))
    return cases


CASES = load_cases()


def divergent_cases() -> dict[str, str]:
    """Cases where the port deliberately differs from pal2nal.v14.pl.

    Maps case name -> reason. For these, the reference lives in expected/
    rather than golden/, so the Perl's behaviour stays on record in
    golden/ and the divergence is visible as a separate, reviewable file.
    """
    path = TESTS / "divergences.tsv"
    if not path.exists():
        return {}
    out = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name, _, reason = line.partition("\t")
        out[name] = reason
    return out


DIVERGENT = divergent_cases()


def _read(path: Path) -> str:
    """A golden is captured bytes, not necessarily valid UTF-8: v14 echoed
    whatever it was given, so a case built on undecodable input has
    undecodable output. Surrogates keep the comparison byte-exact."""
    return path.read_text(encoding="utf-8", errors="surrogateescape")


def read_reference(name: str, ext: str) -> str:
    """The output the port must produce: expected/ when the port is meant
    to differ from the Perl, golden/ otherwise."""
    if name in DIVERGENT:
        path = EXPECTED / f"{name}.{ext}"
        if path.exists():
            return _read(path)
        if (EXPECTED / f"{name}.out").exists():
            return ""      # divergence recorded, this stream is empty
    path = GOLDEN / f"{name}.{ext}"
    return _read(path) if path.exists() else ""


@pytest.mark.parametrize("name,args", CASES, ids=[c[0] for c in CASES])
def test_matches_perl(name: str, args: list[str], monkeypatch: pytest.MonkeyPatch) -> None:
    from pal2nal.cli import main

    # the goldens embed the input paths, so run from the tests directory
    monkeypatch.chdir(TESTS)

    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        status = main(args)

    why = DIVERGENT.get(name)
    note = f" (deliberate divergence: {why})" if why else ""
    expected_status = int(read_reference(name, "status").strip() or 0)
    assert out.getvalue() == read_reference(name, "out"), f"stdout differs for {name}{note}"
    assert err.getvalue() == read_reference(name, "err"), f"stderr differs for {name}{note}"
    assert status == expected_status, f"exit status differs for {name}{note}"


def test_every_divergence_is_documented() -> None:
    """A divergence must name a case that exists and carry a reason."""
    names = {c[0] for c in CASES}
    for name, reason in DIVERGENT.items():
        assert name in names, f"divergences.tsv names an unknown case: {name}"
        assert reason.strip(), f"divergence {name} has no stated reason"
        assert (EXPECTED / f"{name}.out").exists(), f"no expected/ output for {name}"
