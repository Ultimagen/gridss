#!/usr/bin/env python3
# DESCRIPTION : Rematch reads to haplotypes in a given assembly file.

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


MATCH_SCORE = 1
MISMATCH_SCORE = -4
GAP_OPEN_PENALTY = -6
GAP_EXTENSION_PENALTY = -1

matrix = parasail.matrix_create("ACGT", MATCH_SCORE, MISMATCH_SCORE)

def align_parasail_local(seq1, seq2, match_score_only):
    """
    Perform a Smith–Waterman (local) alignment with affine gaps.
    Returns a TraceResult object with .score, .cigar, .end_query, .end_ref, etc.
    """
    # sw_trace_striped_16 does 16‑bit SIMD traces; you can also use sw_trace_scan_16
    if match_score_only:
        res = parasail.sw_striped_16(
            seq1, seq2,
            abs(GAP_OPEN_PENALTY),  # Parasail expects positive open/extend values
            abs(GAP_EXTENSION_PENALTY),
            matrix
        )
        return res
    else:
        res = parasail.sw_trace_striped_16(
            seq1, seq2,
            abs(GAP_OPEN_PENALTY),  # Parasail expects positive open/extend values
            abs(GAP_EXTENSION_PENALTY),
            matrix
        )
        return res

def run_alignment(fa_seq, sequence, match_score_only = False):
    """
    Perform the alignment between two sequences
    @param fa_seq: The reference sequence
    @param sequence: The query sequence
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

        start_pos_adjust, end_pos_adjust = adjust_start_end_positions(t_gap, len(alignment.traceback.ref))

        return alignment.score, start_pos_adjust, end_pos_adjust


def adjust_start_end_positions(t_gap, alignment_length):
    """
    Adjust start and end positions according to sc and the beginning of the alignment
    @param t_gap: The gap at the beginning of the alignment
    @param alignment_length: The length of the alignment
    @return: The adjusted start and end positions
    """
    adjusted_start_pos = 0

    if t_gap > 0:
        adjusted_start_pos = t_gap

    adjusted_end_pos = adjusted_start_pos + alignment_length

    return adjusted_start_pos, adjusted_end_pos


def find_best_haplotype(region_haps, read, reference):
    """
    Find the best haplotype for a given read based on local alignment scores.
    @param region_haps: The haplotypes in the region
    @param read: The read to align
    @param reference: The reference genome
    @return: The best haplotype, its score, start and end points of the alignment, and affected haplotypes
    """
    # in case cigar starts with soft clip, we need to adjust the haplotype start
    affected_haps = set()
    sc_length_start = read.cigar[0][1] if read.cigar[0][0] == pysam.CSOFT_CLIP else 0
    sc_length_end = read.cigar[-1][1] if read.cigar[-1][0] == pysam.CSOFT_CLIP else 0

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
        hap_sc_size_start = hap.cigar[0][1] if hap.cigar[0][0] == pysam.CSOFT_CLIP else 0
        hap_sc_size_end = hap.cigar[-1][1] if hap.cigar[-1][0] == pysam.CSOFT_CLIP else 0

        hap_seq = hap.query_sequence
        hap_start_position = hap.reference_start - hap_sc_size_start
        hap_end_position = hap.reference_end + hap_sc_size_end

        # local alignment
        if((read_start_position <= hap_end_position) and
            (read_end_position >= hap_start_position)):

            score, _, _ = run_alignment(hap_seq,
                                          read_seq,
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
                                      True)
    if ref_score > best_score:
        best_score = ref_score
        best_hap = None
    else:
        # Rerun best haplotype alignment with the best score
        score, start_pos_local, end_pos_local = run_alignment(best_hap_seq,
                                      read_seq,
                                      False)

        best_start_point = start_pos_local
        best_end_point = end_pos_local


    # in case hap direction is reverse, we need to reverse the start and end points
    if best_hap is not None and best_hap.is_reverse:
        best_start_point, best_end_point = max(best_hap.query_length - best_end_point + 1,
                                               0), best_hap.query_length - best_start_point + 1


    return best_hap, best_score, best_start_point, best_end_point, affected_haps

def count_mismatches(md_tag):
    # Remove deletion segments (^...) from MD string
    md_clean = re.sub(r'\^[A-Z]+', '', md_tag)
    # Find all letters that represent mismatches
    mismatches = re.findall(r'[A-Z]', md_clean)
    return len(mismatches)

def rematch_reads_to_haplotypes_in_contig(assembly_path, tumor_crams, germline_crams, reference_path, bed_file_regions, min_sc_indel_size, min_mismatch_count, min_mapq, contig, output_path):
    """
    Rematch reads to haplotypes in a given assembly file.
    @param assembly_path: The input assembly file
    @param tumor_crams: The input tumor CRAM files
    @param germline_crams: The input germline CRAM files
    @param reference_path: The reference genome FASTA file
    @param bed_file_regions: The bed file with the regions to realign
    @param min_sc_indel_size: Minimum size of an indel and soft-clipping in the read to include the read in the assembly
    @param min_mismatch_count: Minimal number of counts to require to include the read in the assembly
    @param min_mapq: Minimum mapping quality to consider a read
    @param contig: The contig to process
    @param output_path: The output assembly file with realigned supporting reads to haplotypes
    """
    logger.info(f"Rematching reads to haplotypes on contig: {contig}")
    reference = pyfaidx.Fasta(reference_path, build_index=False)

    haps_map = dict()
    total_affected_haps = set()
    # align tumor and germline reads to haplotype
    logger.debug("Aligning reads to haplotypes")
    if tumor_crams and germline_crams:
        crams_array = [(cram, 0) for cram in (tumor_crams or [])] + [(cram, 1) for cram in (germline_crams or [])]
    else: # if only one type of CRAM files is provided
        crams_array = [(cram, 0) for cram in (tumor_crams or [])] + [(cram, 0) for cram in (germline_crams or [])]

    min_sc_indel_size_values = [int(x) for x in min_sc_indel_size.split(";")] if min_sc_indel_size else None
    min_mismatch_count_values = [int(x) for x in min_mismatch_count.split(";")] if min_mismatch_count else None
    for cram_file, category in crams_array:
        logger.debug(f"Processing {cram_file} with category {category}")
        # open the CRAM file for fetching the reads
        with pysam.AlignmentFile(cram_file, reference_filename=reference_path) as reads_cram:
            # open assembly file for fetching the haplotypes
            with pysam.AlignmentFile(assembly_path, "rb") as assembly:
                with open(bed_file_regions, "r") as bed:
                    for line in bed:
                        logger.debug(f"Processing line: {line.strip()}")
                        chrom, start, end = line.strip().split()[:3]
                        if chrom != contig:
                            continue
                        start, end = int(start), int(end)
                        # fetch the haplotypes in the region
                        region_haps = [hap for hap in list(assembly.fetch(chrom, start, end))]
                        for read in reads_cram.fetch(chrom, start, end):

                            mapq = read.mapping_quality
                            # The MD tag encodes the positions of mismatches in the alignment compared to the reference
                            if read.has_tag("MD"):
                                md = read.get_tag("MD")
                                mismatch_count = count_mismatches(md)
                            else:
                                mismatch_count = None
                            # Exclude PCR/optical duplicates, low mapping quality reads, and reads with small number of mismatches
                            if (not read.is_duplicate) and (mapq > min_mapq) and \
                                    ((not min_sc_indel_size_values and not min_mismatch_count_values) or
                                    (min_sc_indel_size_values and (any(op in {1, 2, 4} and
                                                                        length > min_sc_indel_size_values[category] for
                                                                        op, length in (read.cigartuples or [])))) or
                                    (min_mismatch_count_values and (
                                            mismatch_count is not None and mismatch_count >= min_mismatch_count_values[category]))):
                                logger.debug(f"Processing read: {read.query_name} with cigartuples {read.cigartuples} ")
                                # Find the best haplotype for the read
                                best_hap, best_score, start_point, end_point, affected_haps = find_best_haplotype(region_haps,
                                                                                                                  read,
                                                                                                                  reference)
                                # If a best haplotype was found, add it to the haps_map
                                if best_hap is not None:
                                    if (best_hap.query_name, best_hap.flag) not in haps_map:
                                        haps_map[(best_hap.query_name, best_hap.flag)] = []
                                    # Append the SupportingRead object to the list
                                    haps_map[(best_hap.query_name, best_hap.flag)].append(
                                        SupportingRead(read.query_name, start_point, end_point, best_score, category))
                                    # update the affected haplotypes
                                    total_affected_haps.update(affected_haps)
                                    logger.debug(f"Found best haplotype: {best_hap.query_name} with score {best_score} and start/end points {start_point}/{end_point}")


    # Write the haplotypes to the output file
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


def parse_args():
    parser = argparse.ArgumentParser(description='Rematch reads to haplotypes')
    parser.add_argument('--assembly', required=True, type=str, help='The input assembly file (BAM format)')
    parser.add_argument('--reference', required=True, type=str, help='The reference genome FASTA file')
    parser.add_argument('--output', required=True, type=str, help='The output assembly file with realigned supporting reads to haplotypes')
    parser.add_argument("--n_jobs", required=False, type=int, default=-1, help="n_jobs of parallel on contigs")
    parser.add_argument("--tumor_crams", required=False, nargs='+', type=str, default=None, help="The input tumor CRAM files")
    parser.add_argument("--germline_crams", required=False, nargs='+', type=str, default=None, help="The input germline CRAM files")
    parser.add_argument("--bed_file_regions", required=True, type=str, help="The bed file with the regions to realign")
    parser.add_argument("--min_mapq", required=True, type=int, help="Minimum mapping quality")
    parser.add_argument("--min_sc_indel_size", required=False, type=str, help="Minimum size of an indel and soft-clipping in the read to include the read in the assembly. ;-separated between samples'")
    parser.add_argument("--min_mismatch_count", required=False, type=str, help="Minimal number of counts to require to include the read in the assembly. ;-separated between samples")


    args = parser.parse_args()
    return args


def run():
    args = parse_args()
    logger.info(f"Running rematch_reads_to_haplotypes_in_contig with args: {args}")

    # Check if the assembly file exists
    if not os.path.exists(args.assembly):
        logger.error(f"Assembly file {args.assembly} does not exist.")
        return

    # Check if the reference file exists
    if not os.path.exists(args.reference):
        logger.error(f"Reference file {args.reference} does not exist.")
        return


    # Check if the bed file exists
    if not os.path.exists(args.bed_file_regions):
        logger.error(f"Bed file {args.bed_file_regions} does not exist.")
        return

    # check if the tumor CRAM files exist
    if args.tumor_crams:
        for cram in args.tumor_crams:
            if not os.path.exists(cram):
                logger.error(f"Tumor CRAM file {cram} does not exist.")
                return

    # check if the germline CRAM files exist
    if args.germline_crams:
        for cram in args.germline_crams:
            if not os.path.exists(cram):
                logger.error(f"Germline CRAM file {cram} does not exist.")
                return

    # Check if the minimum soft-clipping and indel size is valid
    min_sc_indel_size = args.min_sc_indel_size.split(";") if args.min_sc_indel_size else None
    if min_sc_indel_size and (( args.tumor_crams and args.germline_crams and len(min_sc_indel_size) != 2) or
            ((not args.tumor_crams or not args.germline_crams) and len(min_sc_indel_size) != 1) or
            not all(x.isdigit() for x in min_sc_indel_size) or
            not all(int(x) >= 0 for x in min_sc_indel_size)):
        logger.error(f"Minimum soft-clipping and indel size should be a semicolon-separated list of positive integers: {args.min_sc_indel_size}.")
        return

    # Check if the minimum mismatch count is valid
    min_mismatch_count = args.min_mismatch_count.split(";") if args.min_mismatch_count else None
    if min_mismatch_count and ((args.tumor_crams and args.germline_crams and len(min_mismatch_count) != 2) or
            ((not args.tumor_crams or not args.germline_crams) and len(min_mismatch_count) != 1) or
            not all(x.isdigit() for x in min_mismatch_count) or
            not all(int(x) >= 0 for x in min_mismatch_count)):
        logger.error(f"Minimum mismatch count should be a semicolon-separated list of positive integers: {args.min_mismatch_count}.")
        return

    # Check if the minimum mapping quality is valid
    if not isinstance(args.min_mapq, int) or args.min_mapq < 0:
        logger.error(f"Minimum mapping quality should be a positive integer: {args.min_mapq}.")
        return

    with pysam.AlignmentFile(args.assembly, "rc") as assembly_file:
        # Get the list of contig names
        contigs = assembly_file.references
        # Get the list of contig lengths
        contig_lengths = assembly_file.lengths
        large_contigs = [
            contigs[i] for i in range(len(contigs)) if contig_lengths[i] > MIN_CONTIG_LENGTH
        ]

        results = Parallel(n_jobs=args.n_jobs, backend="multiprocessing", max_nbytes=None)(
            delayed(rematch_reads_to_haplotypes_in_contig)(
                args.assembly,
                args.tumor_crams,
                args.germline_crams,
                args.reference,
                args.bed_file_regions,
                args.min_sc_indel_size,
                args.min_mismatch_count,
                args.min_mapq,
                contig,
                f"{args.output}{contig}.bam"
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


if __name__ == "__main__":
    run()