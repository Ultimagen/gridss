import random
import sys
from pathlib import Path

import pysam
import pytest
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import align_long_homopolymers as alh  # noqa: E402

# Same parameters as realign_homopolymers
GLOBAL = alh.create_aligner("global", 1, -4, -6, -1, -5)
LOCAL = alh.create_aligner("local", 1, -4, -6, -1, 0)
REF = "GCTAAAGACAATTACATAACATACACGTCA"
REF_START = 100
SC_LENGTH = 2


def rc(seq):
    return str(Seq(seq).reverse_complement())


def make_read(seq):
    read = pysam.AlignedSegment()
    read.query_sequence = seq
    read.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + i % 40) for i in range(len(seq))))
    return read


def align_and_choose_four_alignments(fa_seq, query):
    """The choice before score-first: all four alignments built, then compared"""
    best = []
    for seq in (query, rc(query)):
        glob = alh.run_alignment(fa_seq, seq, REF_START, SC_LENGTH, GLOBAL)
        local = alh.run_alignment(fa_seq, seq, REF_START, SC_LENGTH, LOCAL)
        if glob[0] >= local[0] + GLOBAL.query_end_gap_score * local[1].count("S"):
            best.append((glob, True))
        else:
            best.append((local, False))
    (orig, orig_is_global), (rev, rev_is_global) = best
    if orig[0] >= rev[0]:
        return 0 if orig_is_global else 1, orig[1], orig[2], orig[4]
    return 2 if rev_is_global else 3, rev[1], rev[2], rev[4]


@pytest.mark.parametrize(
    "fa_seq, query, expected_choice, expected_cigar",
    [
        pytest.param(REF, REF, 0, "30M", id="original_global"),
        pytest.param(REF, REF[:22] + "TTTTTTTT", 1, "22M8S", id="original_local"),
        pytest.param(REF, rc(REF), 2, "30M", id="reverse_global"),
        pytest.param(REF, rc(REF[:22] + "TTTTTTTT"), 3, "22M8S", id="reverse_local"),
        # A reverse local score (10) above the original (9) but the soft clip penalty picks the reverse global (0)
        pytest.param("AGCTTCATGAAGTC", "AGCTTCATGAAGTAAAGGA", 0, "14M5S", id="reverse_bound_loose"),
        pytest.param("", "ACGT", 0, "", id="empty_reference"),
    ],
)
def test_align_and_choose_picks_the_winner(fa_seq, query, expected_choice, expected_cigar):
    read = make_read(query)
    _, cigar, start_pos, r_start, choice, _ = alh.align_and_choose(
        read, query, GLOBAL, LOCAL, fa_seq, query, REF_START, SC_LENGTH, []
    )
    assert (choice, cigar) == (expected_choice, expected_cigar)
    if fa_seq:
        assert (choice, cigar, start_pos, r_start) == align_and_choose_four_alignments(fa_seq, query)


def test_align_and_choose_tie_keeps_the_original_strand():
    palindrome = "ACGTTGCAATGCATTGCAACGT"
    assert palindrome == rc(palindrome)
    read = make_read(palindrome)
    qual = read.qual
    _, cigar, _, _, choice, dels = alh.align_and_choose(
        read, palindrome, GLOBAL, LOCAL, palindrome, palindrome, REF_START, SC_LENGTH, [(3, 5)]
    )
    assert (choice, cigar, dels) == (0, "22M", [(3, 5)])
    assert (read.query_sequence, read.qual) == (palindrome, qual)


def test_align_and_choose_reverse_complements_the_read():
    sequence = rc(REF)[:10] + "A" * 12 + rc(REF)[10:20] + "C" * 11 + rc(REF)[20:]
    edited_sequence, dels, _ = alh.remove_long_homopolymers(sequence)
    read = make_read(sequence)
    qual = read.qual
    _, _, _, _, choice, rev_dels = alh.align_and_choose(
        read, sequence, GLOBAL, LOCAL, REF, edited_sequence, REF_START, SC_LENGTH, dels
    )
    assert choice in (2, 3)
    assert read.query_sequence == rc(sequence)
    assert read.qual == qual[::-1]
    assert dels == [(20, 22), (42, 43)]
    assert rev_dels == [(10, 11), (31, 33)]


def test_align_and_choose_matches_four_alignments_on_random_input():
    rng = random.Random(0)
    for _ in range(2000):
        fa_seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 40)))
        source = rng.choice([fa_seq, rc(fa_seq), fa_seq[: len(fa_seq) // 2] + rc(fa_seq[len(fa_seq) // 2 :])])
        query = list(source)
        for _ in range(rng.randint(0, 4)):
            query.insert(rng.randrange(len(query) + 1), rng.choice("ACGT"))
        query = "".join(query)
        _, cigar, start_pos, r_start, choice, _ = alh.align_and_choose(
            make_read(query), query, GLOBAL, LOCAL, fa_seq, query, REF_START, SC_LENGTH, []
        )
        assert (choice, cigar, start_pos, r_start) == align_and_choose_four_alignments(fa_seq, query), (fa_seq, query)


def test_local_score_bounds_the_global_score():
    rng = random.Random(1)
    for _ in range(5000):
        fa_seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 40)))
        query = "".join(rng.choice([fa_seq[i : i + 5], rng.choice("ACGT") * rng.randint(1, 3)]) for i in range(0, 40, 5))
        query = query or "A"
        assert alh.alignment_score(fa_seq, query, LOCAL) >= alh.alignment_score(fa_seq, query, GLOBAL), (fa_seq, query)


def remove_long_homopolymers_by_loop(sequence, homopolymer_length):
    """The character loop the regex replaced"""
    i, del_length, edited_sequence, start_end_del_tuples = 0, 0, "", []
    while i < len(sequence):
        count = 1
        while i + count < len(sequence) and sequence[i] == sequence[i + count]:
            count += 1
        if count > homopolymer_length:
            edited_sequence += sequence[i] * homopolymer_length
            i += count
            start_end_del_tuples.append((i - (count - homopolymer_length), i))
            del_length += count - homopolymer_length
        else:
            edited_sequence += sequence[i] * count
            i += count
    return edited_sequence, start_end_del_tuples, del_length


def test_remove_long_homopolymers():
    assert alh.remove_long_homopolymers("CA" + "T" * 13 + "GG" + "C" * 10 + "A", 10) == (
        "CA" + "T" * 10 + "GG" + "C" * 10 + "A",
        [(12, 15)],
        3,
    )
    assert alh.remove_long_homopolymers("", 10) == ("", [], 0)


@pytest.mark.parametrize("homopolymer_length", [1, 3, 10])
def test_remove_long_homopolymers_matches_loop(homopolymer_length):
    rng = random.Random(homopolymer_length)
    for _ in range(500):
        sequence = "".join(
            rng.choice("ACGTacgtN") * rng.choice([1, 1, 2, 4, 11, 15]) for _ in range(rng.randint(0, 12))
        )
        assert alh.remove_long_homopolymers(sequence, homopolymer_length) == remove_long_homopolymers_by_loop(
            sequence, homopolymer_length
        )


@pytest.mark.parametrize(
    "cigar, expected",
    [
        ("10S90M", False),
        ("5S20M3I10M2D60M5S", False),
        ("20M3D2I75M", True),
        ("20M3I2I75M", True),
        ("50M50M", True),
        ("10H90M", False),
    ],
)
def test_cigar_may_change(cigar, expected):
    assert alh.cigar_may_change(alh.cigar_string_to_cigartuples(cigar)) == expected


def test_cigar_may_change_false_only_when_adjust_is_a_no_op():
    rng = random.Random(0)
    for _ in range(2000):
        cigar = "".join(f"{rng.randint(1, 20)}{rng.choice('MIDSH')}" for _ in range(rng.randint(1, 6)))
        if not alh.cigar_may_change(alh.cigar_string_to_cigartuples(cigar)):
            assert alh.adjust_and_merge_cigar(cigar) == cigar, cigar


def test_ambiguous_base_pattern_matches_the_old_scan():
    for sequence in ["ACGT", "acgtACGT", "ACGNT", "ACGT.", "", "RYK"]:
        old_skip = not all(c in "atgc" for c in sequence.lower())
        assert (alh.AMBIGUOUS_BASE.search(sequence) is not None) == old_skip
