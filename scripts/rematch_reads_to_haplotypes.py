# DESCRIPTION
#    This script realigns haplotypes in the areas that contain long homopolymer runs
#
from itertools import accumulate
import numpy as np
import pysam
import argparse
import logging
import subprocess
import array
import pyfaidx
from joblib import Parallel, delayed
import os
import parasail
import re


logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__ if __name__ != "__main__" else "rematch_reads_to_haplotypes")

MIN_CONTIG_LENGTH = 100000


class SupportingRead:
    def __init__(self, read_name, start_overlap, end_overlap, score, category):
        self.read_name = read_name
        self.start_overlap = start_overlap
        self.end_overlap = end_overlap
        self.score = score
        self.category = category


match = 1
mismatch = -4
gap_penalty = -6
gap_extension_penalty = -1

matrix = parasail.matrix_create("ACGT", match, mismatch)

def align_parasail_local(seq1, seq2, match_score_only):
    """
    Perform a Smith–Waterman (local) alignment with affine gaps.
    Returns a TraceResult object with .score, .cigar, .end_query, .end_ref, etc.
    """
    # sw_trace_striped_16 does 16‑bit SIMD traces; you can also use sw_trace_scan_16
    if match_score_only:
        res = parasail.sw_striped_16(
            seq1, seq2,
            abs(gap_penalty),  # Parasail expects positive open/extend values
            abs(gap_extension_penalty),
            matrix
        )
        return res
    else:
        res = parasail.sw_trace_striped_16(
            seq1, seq2,
            abs(gap_penalty),  # Parasail expects positive open/extend values
            abs(gap_extension_penalty),
            matrix
        )
        return res

def run_alignment(fa_seq, sequence, start_pos, sc_length, hap_cigar, match_score_only = False):
    """
    Perform the alignment between two sequences
    @param fa_seq: The reference sequence
    @param sequence: The query sequence
    @param start_pos: The start position of the read
    @param sc_length: The length of soft clipping
    @param hap_cigar: The CIGAR string of the haplotype
    @param match_score_only: If True, return only the match score, for running faster
    @return: The alignment score, CIGAR string, start position, and query start position
    """

    if len(fa_seq) == 0 or len(sequence) == 0:
        return 0, 0, 0

    if match_score_only:
        # for saving time, get score only
        return align_parasail_local(fa_seq, sequence, match_score_only).score, 0, 0
    else:
        # calculate start and end as well
        alignment = align_parasail_local(fa_seq, sequence, match_score_only)
        cigar = alignment.cigar.decode.decode('ascii')
        parsed_cigar = re.match(r'^(\d+)([IS])', cigar)
        t_gap = int(parsed_cigar.group(1)) if parsed_cigar else 0
        t_gap += alignment.cigar.beg_query

        start_pos_adjust, end_pos_adjust = adjust_start_end_positions(start_pos - sc_length, t_gap, len(alignment.traceback.ref), hap_cigar)

        return alignment.score, start_pos_adjust, end_pos_adjust


def adjust_start_end_positions(start_pos, t_gap, alignmnet_length, hap_cigar_tuples):
    """
        Adjust start and end positions according to sc and the beginning of the alignment
    """
    adjusted_start_pos = 0

    if t_gap > 0:
        adjusted_start_pos = t_gap

    adjusted_end_pos = adjusted_start_pos + alignmnet_length

    return start_pos + adjusted_start_pos, start_pos + adjusted_end_pos


def find_best_haplotype(region_haps, read, reference):
    # in case cigar starts with soft clip, we need to adjust the haplotype start
    affected_haps = set()
    sc_length_start = read.cigar[0][1] if read.cigar[0][0] == 4 else 0
    sc_length_end = read.cigar[-1][1] if read.cigar[-1][0] == 4 else 0

    # find the best alignment
    best_score = -np.inf
    best_hap = None
    best_start_point = None
    best_end_point = None
    read_seq = read.query_sequence
    best_hap_start_position = None
    best_hap_seq = None

    read_start_position = read.reference_start - sc_length_start
    read_end_position = read.reference_end + sc_length_end
    # extract the region of the assembly that the read overlaps
    for hap in region_haps:
        # only in case we overlap the breakpoint
        # read.reference_start is the breakpoint position
        hap_sc_size_start = hap.cigar[0][1] if hap.cigar[0][0] == 4 else 0
        hap_sc_size_end = hap.cigar[-1][1] if hap.cigar[-1][0] == 4 else 0

        hap_seq = hap.query_sequence
        hap_start_position = hap.reference_start - hap_sc_size_start
        hap_end_position = hap.reference_end + hap_sc_size_end

        # local alignment
        if((read_start_position <= hap_end_position) and
            (read_end_position >= hap_start_position) and
            ((sc_length_start == 0 and sc_length_end == 0) or (sc_length_start > 0 and read.reference_start >= hap_start_position) or
            (sc_length_end > 0 and read.reference_end <= hap_end_position))):

            score, _, _ = run_alignment(hap_seq,
                                          read_seq,
                                          hap_start_position,
                                          0,
                                          hap.cigartuples,
                                          True)
            if score > best_score:
                best_score = score
                best_hap = hap
                best_hap_seq = hap_seq
                best_hap_start_position = hap_start_position

                affected_haps.add((hap.query_name, hap.flag))

    # check reference sequence as well
    del_length = sum([length for op, length in read.cigartuples if op == 2])

    ref_seq = reference[read.reference_name][
             max(read.reference_start - sc_length_start, 0):
             min(read.reference_end + sc_length_end + del_length, len(reference[read.reference_name]))].seq.upper()
    ref_score, _, _ = run_alignment(ref_seq,
                                    read_seq,
                                    read_start_position,
                                    0,
                                    [],
                                      True)
    if ref_score > best_score:
        best_score = ref_score
        best_hap = None
    else:
        # Rerun best haplotype alignment with the best score
        score, start_pos_local, end_pos_local = run_alignment(best_hap_seq,
                                      read_seq,
                                      best_hap_start_position,
                                      0,
                                      best_hap.cigartuples,
                                       False)

        best_start_point = start_pos_local - best_hap_start_position
        best_end_point = end_pos_local - best_hap_start_position


    # in case hap direction is reverse, we need to reverse the start and end points
    if best_hap is not None and best_hap.is_reverse:
        best_start_point, best_end_point = max(best_hap.query_length - best_end_point + 1,
                                               0), best_hap.query_length - best_start_point + 1


    return best_hap, best_score, best_start_point, best_end_point, affected_haps

def rematch_homopolymere(assembly_path, tumor_crams, germline_crams, reference_path, bed_file_regions, contig, output_path):

    logger.info(f"Rematching reads to haplotypes on contig: {contig}")
    reference = pyfaidx.Fasta(reference_path, build_index=False)

    haps_map = dict()
    total_affected_haps = set()
    # align tumor and germline reads to haplotype
    logger.debug("Aligning reads to haplotypes")
    if tumor_crams:
        crams_array = [(cram, 0) for cram in (tumor_crams or [])] + [(cram, 1) for cram in (germline_crams or [])]
    else:
        crams_array = [(cram, 0) for cram in (germline_crams or [])]
    for cram_file, category in crams_array:
        logger.debug(f"Processing {cram_file} with category {category}")
        with pysam.AlignmentFile(cram_file, reference_filename=reference_path) as reads_cram:
            with open(bed_file_regions, "r") as bed:
                for line in bed:
                    logger.debug(f"Processing line: {line.strip()}")
                    chrom, start, end = line.strip().split()[:3]
                    if chrom != contig:
                        continue
                    start, end = int(start), int(end)
                    # fetch the haplotypes in the region
                    with pysam.AlignmentFile(assembly_path, "rb") as assembly:
                        region_haps = [
                            hap for hap in list(assembly.fetch(chrom, start, end))
                            if any(op in {1, 2, 4} and length > 20 for op, length in (hap.cigartuples or []))
                        ]
                        for read in reads_cram.fetch(chrom, start, end):
                            # Exclude PCR/optical duplicates, keep only tumor reads with long insertions, deletions, or soft-clips
                            if (not read.is_duplicate) and (category == 1 or any(
                                    op in {1, 2, 4} and length > 20 for op, length in (read.cigartuples or []))):
                                logger.debug(f"Processing read: {read.query_name} with cigartuples {read.cigartuples} ")
                                best_hap, best_score, start_point, end_point, affected_haps = find_best_haplotype(region_haps,
                                                                                                                  read,
                                                                                                                  reference)
                                if best_hap is not None:
                                    if (best_hap.query_name, best_hap.flag) not in haps_map:
                                        haps_map[(best_hap.query_name, best_hap.flag)] = []
                                    # Append the SupportingRead object to the list
                                    haps_map[(best_hap.query_name, best_hap.flag)].append(
                                        SupportingRead(read.query_name, start_point, end_point, best_score, category))
                                    # update the affected haplotypes
                                    total_affected_haps.update(affected_haps)
                                    logger.debug(f"Found best haplotype: {best_hap.query_name} with score {best_score} and start/end points {start_point}/{end_point}")


    # write the haplotypes to the output file
    logger.info(f"Writing the haplotypes to the output file contig: {contig}")
    with pysam.AlignmentFile(assembly_path) as assembly:
        with pysam.AlignmentFile(output_path + "_unsorted.bam", "wb", template=assembly) as output:
            # write each haplotype as a row in the output file
            # supporting reads are stored in the ef tag
            for hap in assembly.fetch(contig):
                if (hap.query_name, hap.flag) in haps_map:
                    # in case the haplotype is affected and has reads supporting it
                    supporting_reads = haps_map[(hap.query_name, hap.flag)]
                    hap.set_tag("ef", " ".join([read.read_name for read in supporting_reads]))
                    hap.set_tag("ez", " ".join([read.read_name for read in supporting_reads]))
                    hap.set_tag("eq", array.array("f", [read.score for read in supporting_reads]))
                    hap.set_tag("os", array.array("i", [read.start_overlap for read in supporting_reads]))
                    hap.set_tag("oe", array.array("i", [read.end_overlap for read in supporting_reads]))
                    hap.set_tag("ec", array.array("i", [read.category for read in supporting_reads]))
                    hap.set_tag("et", array.array("b", [0 for read in supporting_reads]))
                    output.write(hap)
                elif (hap.query_name, hap.flag) in total_affected_haps:
                    # in case the haplotype is affected but no reads are supporting it
                    continue
                else:
                    # in case the haplotype is not affected
                    output.write(hap)

    # index output file
    if os.path.exists(f"{output_path}_unsorted.bam"):
        subprocess.check_call(f"samtools sort {output_path}_unsorted.bam -o {output_path}", shell=True)
        subprocess.check_call(f"samtools index {output_path}", shell=True)
        # remove unsoreted file
        subprocess.check_call(f"rm {output_path}_unsorted.bam", shell=True)



parser = argparse.ArgumentParser(description='Rematch reads to haplotypes')
parser.add_argument('--assembly', required=True, help='The input assembly file')
parser.add_argument('--reference', required=True, help='The reference genome FASTA file')
parser.add_argument('--output', required=True, help='The output assembly file with realigned supporting reads to haplotypes')
parser.add_argument("--n_jobs", help="n_jobs of parallel on contigs", type=int, default=-1)
parser.add_argument("--tumor_crams", help="The input tumor CRAM files", nargs='+', required=False)
parser.add_argument("--germline_crams", help="The input germline CRAM files", nargs='+', required=True)
parser.add_argument("--bed_file_regions", help="The bed file with the regions to realign", required=True)
args = parser.parse_args()

# rematch_homopolymere(args.assembly, args.tumor_crams, args.germline_crams, args.reference, args.bed_file_regions, "chr1", args.output)
with pysam.AlignmentFile(args.assembly, "rc") as assembly_file:
    # Get the list of contig names
    contigs = assembly_file.references
    # Get the list of contig lengths
    contig_lengths = assembly_file.lengths
    large_contigs = [
        contigs[i] for i in range(len(contigs)) if contig_lengths[i] > MIN_CONTIG_LENGTH
    ]

    results = Parallel(n_jobs=args.n_jobs, backend="multiprocessing", max_nbytes=None)(
        delayed(rematch_homopolymere)(
            args.assembly, args.tumor_crams, args.germline_crams, args.reference, args.bed_file_regions, contig, f"{args.output}{contig}.bam"
        )
        for contig in large_contigs
    )

    # merge the contig files together
    with pysam.AlignmentFile(args.output, mode='wb', header=assembly_file.header) as output:
        for contig in contigs:
            if contig in large_contigs:
                with pysam.AlignmentFile(f"{args.output}{contig}.bam") as contig_file:
                    for read in contig_file:
                        output.write(read)
            else:
                # copy the reads from the original file in case the contig is short
                for read in assembly_file.fetch(contig):
                        output.write(read)

        # remove the contig files
        for contig in large_contigs:
            os.remove(f"{args.output}{contig}.bam")
            os.remove(f"{args.output}{contig}.bam.bai")

pysam.index(args.output)
