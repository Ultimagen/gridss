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


def test_cigar_may_change_on_tuples_the_string_parser_cannot_produce():
    assert alh.cigar_may_change([])
    assert alh.cigar_may_change([(0, 10), (9, 5), (0, 10)])


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


def write_test_bam(tmp_path, rng):
    contigs = [("chrA", 150_000), ("chrS1", 2_000), ("chrS2", 3_000), ("chrB", 120_000), ("chrS3", 1_000)]
    reference = {name: "".join(rng.choice("ACGT") for _ in range(length)) for name, length in contigs}
    for name in reference:
        pos = rng.randrange(1_000, len(reference[name]) - 100) if len(reference[name]) > 2_000 else 500
        reference[name] = reference[name][:pos] + "T" * 14 + reference[name][pos + 14 :]
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in reference.items()))
    pysam.faidx(str(fasta))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": n, "LN": ln} for n, ln in contigs]}
    reads = []
    for tid, (name, length) in enumerate(contigs):
        for i in range(60):
            start = rng.randrange(0, length - 200)
            seq = reference[name][start : start + 120]
            read = pysam.AlignedSegment()
            read.query_name = f"{name}_{i}"
            read.reference_id, read.reference_start, read.mapping_quality = tid, start, 60
            clip = rng.choice([0, 12, 20])
            seq = "".join(rng.choice("ACGT") for _ in range(clip)) + seq[clip:]
            if rng.random() < 0.5:
                seq = seq[:60] + "A" * 13 + seq[73:]
            read.query_sequence = seq
            read.cigarstring = f"{clip}S{120 - clip}M" if clip else "120M"
            read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            reads.append(read)
    bam = tmp_path / "in.bam"
    with pysam.AlignmentFile(str(tmp_path / "unsorted.bam"), "wb", header=header) as out:
        for read in reads:
            out.write(read)
    pysam.sort("-o", str(bam), str(tmp_path / "unsorted.bam"))
    pysam.index(str(bam))
    return bam, fasta, contigs


def serial_merge(bam, fasta, contigs, tmp_path):
    """The merge before the parallel sort: each large contig realigned and sorted, short contigs copied"""
    records = []
    with pysam.AlignmentFile(str(bam)) as source:
        for name, length in contigs:
            if length > 100_000:
                alh.realign_homopolymers(str(bam), str(tmp_path / f"s_{name}.bam"), str(fasta), 10, name, 10)
                pysam.sort("-o", str(tmp_path / f"s_{name}_sorted.bam"), str(tmp_path / f"s_{name}.bam"))
                with pysam.AlignmentFile(str(tmp_path / f"s_{name}_sorted.bam")) as part:
                    records.extend(read.to_string() for read in part)
            else:
                records.extend(read.to_string() for read in source.fetch(name))
    return records


def test_main_matches_serial_merge(tmp_path, monkeypatch):
    bam, fasta, contigs = write_test_bam(tmp_path, random.Random(3))
    output = tmp_path / "out.bam"
    monkeypatch.setattr(
        sys, "argv", ["alh", "--input", str(bam), "--output", str(output), "--reference", str(fasta), "--n_jobs", "2"]
    )
    alh.main()
    with pysam.AlignmentFile(str(output)) as result, pysam.AlignmentFile(str(bam)) as source:
        assert str(result.header) == str(source.header)
        records = [read.to_string() for read in result]
    assert records == serial_merge(bam, fasta, contigs, tmp_path)
    assert sorted(p.name for p in tmp_path.glob("out.bam*")) == ["out.bam", "out.bam.bai"]
