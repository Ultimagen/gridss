import sys
from pathlib import Path

import pysam
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import choose_best_haplotype_realignment as cbh  # noqa: E402

HEADER = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": 10_000}, {"SN": "chr2", "LN": 10_000}]}


def write_source(path, alignments):
    """alignments: (name, cigar, mapq, flag) per record, in queryname order"""
    with pysam.AlignmentFile(str(path), "wb", header=HEADER) as out:
        for name, cigar, mapq, flag in alignments:
            read = pysam.AlignedSegment()
            read.query_name, read.flag, read.reference_id, read.reference_start = name, flag, 0, 100
            read.mapping_quality, read.cigarstring = mapq, cigar
            read.query_sequence = "A" * 100
            out.write(read)


@pytest.mark.parametrize("uncompressed", [False, True])
def test_run_picks_the_best_source(tmp_path, monkeypatch, uncompressed):
    # r1: source 0 spans more but its supplementary is poor, so source 1 wins; r2: tie, source 0 kept; r3: source 0 spans more
    write_source(tmp_path / "s0.bam", [
        ("r1", "50M50S", 60, 0), ("r1", "50S50M", 10, 2048), ("r2", "100M", 60, 0), ("r3", "90M10S", 60, 0)])
    write_source(tmp_path / "s1.bam", [("r1", "60M40S", 60, 0), ("r2", "100M", 60, 0), ("r3", "80M20S", 60, 0)])
    output = tmp_path / "out.bam"
    argv = ["cbh", "--alignment_source", str(tmp_path / "s0.bam"), "--alignment_source", str(tmp_path / "s1.bam"),
            "--output", str(output), "--min_mapping_quality", "60"]
    monkeypatch.setattr(sys, "argv", argv + (["--uncompressed"] if uncompressed else []))
    cbh.run()
    with pysam.AlignmentFile(str(output)) as result:
        assert [(r.query_name, r.cigarstring) for r in result] == [("r1", "60M40S"), ("r2", "100M"), ("r3", "90M10S")]
    with open(output, "rb") as fh:
        # BGZF header, then deflate block type: stored (level 0) for uncompressed output
        block = fh.read(19)
    assert (block[18] & 0b110 == 0) == uncompressed
