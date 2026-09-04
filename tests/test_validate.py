"""The input gates in pal2nal/validate.py.

The property that matters most is negative: a message about rejected input
must not contain the input. v14's equivalent failure echoed the sequence,
and under -html that echo was the page.
"""

from __future__ import annotations

import pytest

from pal2nal.validate import (
    InputError,
    check_alignment_widths,
    check_nucleotide_alphabet,
    check_text,
)


def message_of(text: str, source: str = "input") -> str:
    with pytest.raises(InputError) as exc:
        check_text(text, source)
    return exc.value.message


def test_plain_ascii_passes() -> None:
    check_text(">s1\nMAKQL-RT\n\t\r\n", "input")


@pytest.mark.parametrize(
    "offender",
    [
        "–",          # en dash, from a word processor
        "﻿",          # UTF-8 byte order mark
        "é",          # latin-1 letter
        "\udce9",          # an undecodable byte, via surrogateescape
        "\u200B",     # zero-width space: invisible, and rejected all the same
    ],
)
def test_non_ascii_is_refused_without_echoing_it(offender: str) -> None:
    message = message_of(f">s1\nMAK{offender}QL\n")
    assert "non-ASCII" in message
    assert offender not in message
    assert "MAK" not in message, "no input text may appear in the message"


def test_non_ascii_message_counts_and_locates() -> None:
    message = message_of("ATG\nATéGé\n")
    assert "contains 2 non-ASCII characters" in message
    assert "line 2, column 3" in message


@pytest.mark.parametrize("newline", ["\n", "\r\n", "\r"])
def test_position_is_counted_for_every_line_ending(newline: str) -> None:
    message = message_of(f"ATG{newline}ATG{newline}ATéG")
    assert "line 3, column 3" in message


def test_tab_cr_and_lf_are_not_control_characters() -> None:
    check_text("s1\tMAKQL\r\nsN\tMAKQL\r", "input")


@pytest.mark.parametrize("code", [0x00, 0x07, 0x0C, 0x1B, 0x7F])
def test_control_characters_are_refused_by_code_point(code: int) -> None:
    message = message_of(f">s1\nMAK{chr(code)}QL\n")
    assert "control character" in message
    assert f"0x{code:02X}" in message
    assert chr(code) not in message


def test_non_ascii_is_reported_before_a_control_character() -> None:
    """Both are refusals; the encoding one is the more useful diagnosis."""
    assert "non-ASCII" in message_of("MAK\x00QLé")


def test_the_source_is_named_verbatim() -> None:
    assert "The protein alignment" in message_of("é", "The protein alignment")


def test_iupac_and_lower_case_dna_pass() -> None:
    check_nucleotide_alphabet({"s1": "ACGTUacgtu", "s2": "RYSWKMBDHVN", "s3": ""})


def test_stray_letters_are_named_with_their_position() -> None:
    with pytest.raises(InputError) as exc:
        check_nucleotide_alphabet({"BC070280": "ACGQTLQ"})
    message = exc.value.message
    assert "BC070280" in message
    assert "'Q' (x2, first at nucleotide 4)" in message
    assert "'L' (x1, first at nucleotide 6)" in message


def test_the_report_is_bounded() -> None:
    """A wholly corrupt file must produce a report, not a second data dump."""
    with pytest.raises(InputError) as exc:
        check_nucleotide_alphabet({f"s{i}": "QQQ" for i in range(20)})
    message = exc.value.message
    assert message.count("first at nucleotide") == 5
    assert "and 15 more sequence(s)" in message


# --------------------------------------------------------------------------
# the alignment must be a rectangle
# --------------------------------------------------------------------------


def test_equal_width_rows_pass() -> None:
    check_alignment_widths(["a", "b", "c"], ["MAKQ", "M-KQ", "MAK."])


@pytest.mark.parametrize("rows", [[], ["MAKQ"]])
def test_nothing_to_compare_passes(rows: list[str]) -> None:
    check_alignment_widths([f"s{i}" for i in range(len(rows))], rows)


def test_a_ragged_row_is_refused() -> None:
    with pytest.raises(InputError) as e:
        check_alignment_widths(["query", "ref"], ["AX1-GTV", "AIHGTV"])
    assert "query: 7" in e.value.message
    assert "ref: 6" in e.value.message


def test_row_order_does_not_change_the_verdict() -> None:
    """v14's answer depended on it: the first row set the alignment length."""
    ids, rows = ["query", "ref"], ["AX1-GTV", "AIHGTV"]
    with pytest.raises(InputError) as first:
        check_alignment_widths(ids, rows)
    with pytest.raises(InputError) as second:
        check_alignment_widths(ids[::-1], rows[::-1])
    assert "query: 7" in second.value.message and "ref: 6" in second.value.message
    assert first.value.message.count("\n") == second.value.message.count("\n")


def test_a_majority_width_is_named_and_only_the_odd_rows_listed() -> None:
    with pytest.raises(InputError) as e:
        check_alignment_widths(["a", "b", "c", "d"], ["MVLSQ", "MVLSQ", "MVLS", "MVLSQ"])
    message = e.value.message
    assert "most have 5 columns" in message
    assert "c: 4" in message
    for absent in ("a: 5", "b: 5", "d: 5"):
        assert absent not in message


def test_no_majority_is_not_claimed_as_one() -> None:
    """Two rows of different lengths give no grounds for calling either intended."""
    with pytest.raises(InputError) as e:
        check_alignment_widths(["a", "b"], ["MVLSQ", "MVLS"])
    assert "most have" not in e.value.message


def test_the_report_is_bounded_and_never_echoes_the_rows() -> None:
    rows = ["MVLSQAD"] * 10 + ["MVLSQADK"] * 8  # majority 7, eight odd rows
    ids = [f"s{i}" for i in range(len(rows))]
    with pytest.raises(InputError) as e:
        check_alignment_widths(ids, rows)
    message = e.value.message
    assert "most have 7 columns" in message
    assert "... and 3 more sequence(s)" in message
    assert "MVLSQAD" not in message
    assert len(message.splitlines()) < 12
