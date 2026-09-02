"""Record the port's output for the cases listed in divergences.tsv.

Run this only after confirming by hand that each difference from golden/ is
explained by an agreed change in PORTING.md. It refuses to touch any case
that is not registered as a divergence.
"""

from __future__ import annotations

import io
import os
import shlex
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

TESTS = Path(__file__).parent
sys.path.insert(0, str(TESTS.parent))

from pal2nal.cli import main  # noqa: E402


def main_() -> None:
    cases = dict(
        line.split("\t", 1) for line in (TESTS / "cases.tsv").read_text().splitlines() if line
    )
    divergent = [
        line.split("\t")[0]
        for line in (TESTS / "divergences.tsv").read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    (TESTS / "expected").mkdir(exist_ok=True)
    os.chdir(TESTS)
    for name in divergent:
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            status = main(shlex.split(cases[name]))
        # surrogateescape: the port's own output is ASCII, but the writer
        # must not choke if a case ever records otherwise
        write = {"encoding": "utf-8", "errors": "surrogateescape"}
        (TESTS / "expected" / f"{name}.out").write_text(out.getvalue(), **write)
        (TESTS / "expected" / f"{name}.status").write_text(f"{status}\n", **write)
        target = TESTS / "expected" / f"{name}.err"
        if err.getvalue():
            target.write_text(err.getvalue(), **write)
        elif target.exists():
            target.unlink()
        print(f"  recorded {name} (status {status})")


if __name__ == "__main__":
    main_()
