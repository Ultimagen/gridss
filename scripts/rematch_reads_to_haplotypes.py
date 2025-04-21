# DESCRIPTION
#    This script realigns haplotypes in the areas that contain long homopolymer runs,
#    where UG data introduces false variation due to the limit on calling homopolymer length
from itertools import accumulate

import numpy as np
import pysam
import argparse
import logging
import subprocess
from Bio import Align
import array
import pyfaidx
from joblib import Parallel, delayed
import os
import parasail


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

def create_aligner(mode, match, mismatch, gap_penalty, gap_extension_penalty, sc_penalty):
    """
    Create a new PairwiseAligner object with the specified parameters
    """
    # Create a new PairwiseAligner object
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    # Set gap scoring
    aligner.open_gap_score = gap_penalty
    aligner.extend_gap_score = gap_extension_penalty

    # Set end gap scores to 0.0 to not penalize them
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = sc_penalty

    return aligner

def run_alignment(fa_seq, sequence, start_pos, sc_length, hap_cigar, aligner, match_score_only = False):
    """
    Perform the alignment between two sequences
    @param fa_seq: The reference sequence
    @param sequence: The query sequence
    @param start_pos: The start position of the read
    @param sc_length: The length of soft clipping
    @param aligner: The alignment object
    @param hap_cigar: The CIGAR string of the haplotype
    @param match_score_only: If True, return only the match score, for running faster
    @return: The alignment score, CIGAR string, start position, and query start position
    """
    # Perform the alignment between two sequences
    if len(fa_seq) == 0 or len(sequence) == 0:
        return 0, 0, 0

    if match_score_only:
        #return next(aligner.align(fa_seq, sequence)).score, 0, 0
        return align_parasail_local(fa_seq, sequence, match_score_only).score, 0, 0
    else:
        # Print alignment's score and the alignment itself
        #alignment_orig = next(aligner.align(fa_seq, sequence))
        alignment = align_parasail_local(fa_seq, sequence, match_score_only)

        start_pos_adjust, end_pos_adjust = adjust_start_end_positions(start_pos - sc_length, alignment.cigar.beg_query, len(alignment.traceback.ref), hap_cigar)
        return alignment.score, start_pos_adjust, end_pos_adjust
    return 0, 0, 0

def adjust_start_end_positions(start_pos, aligned, alignmnet_length, hap_cigar_tuples):
    """
        Convert the aligned segments in biopython format to a CIGAR farmat.
    """
    adjusted_start_pos = 0  # This will store the adjusted start position based on initial target insertions

    if t_gap > 0:
        adjusted_start_pos = t_gap

    adjusted_end_pos = adjusted_start_pos + alignmnet_length
    # adjust end position by hap_cigar_tuples
    accumulate_length = 0
    for op, length in hap_cigar_tuples:


        if op == 1: # insertion
            if accumulate_length <= adjusted_start_pos:
                adjusted_start_pos -= min(length, adjusted_start_pos - accumulate_length)
            if accumulate_length <= adjusted_end_pos:
                adjusted_end_pos -= min(length, adjusted_end_pos - accumulate_length)
            else:
                break

        elif op == 2: # deletion
            if accumulate_length <= adjusted_start_pos:
                adjusted_start_pos += min(length, adjusted_start_pos - accumulate_length)
            if accumulate_length <= adjusted_end_pos:
                adjusted_end_pos += min(length, adjusted_end_pos - accumulate_length)
            else:
                break

        accumulate_length += length

    return start_pos + adjusted_start_pos, start_pos + adjusted_end_pos

def find_best_haplotype(region_haps, read, local_aligner, reference):
    # in case cigar starts with soft clip, we need to adjust the haplotype start
    affected_haps = set()
    sc_size_start = 0
    sc_size_end = 0
    if read.cigar[0][0] == 4:
        read_start = read.reference_start - read.cigar[0][1]
        sc_size_start = read.cigar[0][1]
    else:
        read_start = read.reference_start

    # in case cigar ends with soft clip, we need to adjust the haplotype end
    if read.cigar[-1][0] == 4:
        read_end = read.reference_end + read.cigar[-1][1]
        sc_size_end = read.cigar[-1][1]
    else:
        read_end = read.reference_end

    # find the best alignment
    best_score = -np.inf
    best_hap = None
    best_start_point = None
    best_end_point = None
    read_seq = read.query_sequence
    best_hap_start_position = None
    best_hap_seq = None
    # extract the region of the assembly that the read overlaps
    for hap in region_haps:
        # only in case we overlap the breakpoint
        # read.reference_start is the breakpoint position
        hap_sc_size_start = 0
        hap_sc_size_end = 0
        if hap.cigar[0][0] == 4:
            hap_sc_size_start = hap.cigar[0][1]

        if hap.cigar[-1][0] == 4:
            hap_sc_size_end = hap.cigar[-1][1]

        hap_seq = hap.query_sequence
        hap_start_position = hap.reference_start - hap_sc_size_start
        hap_end_position = hap.reference_end + hap_sc_size_end
        read_start_position = read.reference_start - sc_size_start
        read_end_position = read.reference_end + sc_size_end
        # local alignment
        if((read_start_position <= hap_end_position) and
            (read_end_position >= hap_start_position) and
            ((sc_size_start== 0 and sc_size_end == 0) or (sc_size_start > 0 and read.reference_start >= hap_start_position) or
            (sc_size_end > 0 and read.reference_end <= hap_end_position))):

            score, _, _ = run_alignment(hap_seq,
                                          read_seq,
                                          hap_start_position,
                                          0,
                                          hap.cigartuples,
                                          local_aligner,
                                          True)
            if score > best_score:
                best_score = score
                best_hap = hap
                best_hap_seq = hap_seq
                best_hap_start_position = hap_start_position

                affected_haps.add((hap.query_name, hap.flag))

    # check reference sequence as well
    # sift clip length from the start and end of the read
    sc_length_start = read.cigar[0][1] if read.cigar[0][0] == 4 else 0
    sc_length_end = read.cigar[-1][1] if read.cigar[-1][0] == 4 else 0
    del_length = sum([length for op, length in read.cigartuples if op == 2])

    ref_seq = reference[read.reference_name][
             max(read.reference_start - sc_length_start, 0):
             min(read.reference_end + sc_length_end + del_length, len(reference[read.reference_name]))].seq.upper()
    ref_score, start_pos_local, end_pos_local = run_alignment(ref_seq,
                                                            read_seq,
                                                            read_start,
                                                            sc_size_start,
                                                            [],
                                                            local_aligner,
                                                              False)
    if ref_score > best_score:
        best_score = ref_score
        best_hap = None
        best_start_point = start_pos_local - read_start
        best_end_point = end_pos_local - read_start
    else:
        # Rerun best haplotype alignment with the best score
        score, start_pos_local, end_pos_local = run_alignment(best_hap_seq,
                                      read_seq,
                                      best_hap_start_position,
                                      0,
                                      best_hap.cigartuples,
                                      local_aligner,
                                       False)

        best_hap = hap
        best_start_point = start_pos_local - read_start
        best_end_point = end_pos_local - read_start


    # in case hap direction is reverse, we need to reverse the start and end points
    if best_hap is not None and best_hap.is_reverse:
        best_start_point, best_end_point = max(best_hap.query_length - best_end_point + 1, 0), best_hap.query_length - best_start_point + 1

    return best_hap, best_score, best_start_point, best_end_point, affected_haps

def rematch_homopolymere(assembly_path, tumor_crams, germline_crams, reference_path, bed_file_regions, contig, output):

    logger.info(f"Rematching reads to haplotypes on contig: {contig}")
    local_aligner = create_aligner('local', match, mismatch, gap_penalty, gap_extension_penalty, 0)
    reference = pyfaidx.Fasta(reference_path, build_index=False)

    haps_map = dict()
    total_affected_haps = set()
    # align tumor and germline reads to haplotype
    logger.debug("Aligning reads to haplotypes")
    for cram_file, category in [(cram, 0) for cram in tumor_crams] + [(cram, 1) for cram in germline_crams]:
        logger.debug(f"Processing {cram_file} with category {category}")
        with pysam.AlignmentFile(cram_file) as reads_cram:
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
                                                                                                                  local_aligner,
                                                                                                                  reference)
                                if best_hap is not None:
                                    if best_hap not in haps_map:
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
        with pysam.AlignmentFile(output + "_unsorted.bam", "wb", template=assembly) as output:
            # write each haplotype as a row in the output file
            # supporting reads are stored in the ef tag
            for hap in assembly.fetch(contig):
                if hap in haps_map:
                    # in case the haplotype is affected and has reads supporting it
                    supporting_reads = haps_map[(hap.query_name, hap.flag)]
                    hap.set_tag("ef", " ".join([read.read_name for read in supporting_reads]))
                    hap.set_tag("ez", " ".join([read.read_name for read in supporting_reads]))
                    hap.set_tag("eq", array.array("f", [read.score for read in supporting_reads]))
                    hap.set_tag("os", array.array("i", [read.start_overlap for read in supporting_reads]))
                    hap.set_tag("oe", array.array("i", [read.end_overlap for read in supporting_reads]))
                    hap.set_tag("ec", array.array("i", [read.category for read in supporting_reads]))
                    hap.set_tag("et",
                                array.array("b", [0 for read in supporting_reads]))  # ?? Not sure about what is et
                    output.write(hap)
                elif (hap.query_name, hap.flag) in total_affected_haps:
                    # in case the haplotype is affected but no reads are supporting it
                    hap.set_tag("ef", "")
                    hap.set_tag("ez", "")
                    hap.set_tag("eq", array.array("f", []))
                    hap.set_tag("os", array.array("i", []))
                    hap.set_tag("oe", array.array("i", []))
                    hap.set_tag("ec", array.array("i", []))
                    hap.set_tag("et", array.array("b", []))
                    output.write(hap)
                else:
                    # in case the haplotype is not affected
                    output.write(hap)

    # index output file
    if os.path.exists(f"{output}_unsorted.bam"):
        subprocess.check_call(f"samtools sort {output}_unsorted.bam -o {output}", shell=True)
        subprocess.check_call(f"samtools index {output}", shell=True)
        # remove unsoreted file
        subprocess.check_call(f"rm {output}_unsorted.bam", shell=True)



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
            args.assembly, args.tumor_crams, args.germline_crams, args.reference, args.bed_file_regions, contig, f"{args.output}{contig}"
        )
        for contig in large_contigs
    )

    # merge the contig files together
    with pysam.AlignmentFile(args.output, mode='wb', header=assembly_file.header) as output:
        for contig in contigs:
            if contig in large_contigs:
                with pysam.AlignmentFile(f"{args.output}{contig}_sorted.bam") as contig_file:
                    for read in contig_file:
                        output.write(read)
            else:
                # copy the reads from the original file in case the contig is short
                for read in assembly_file.fetch(contig):
                        output.write(read)

        # remove the contig files
        for contig in large_contigs:
            os.remove(f"{args.output}{contig}_sorted.bam")
            os.remove(f"{args.output}{contig}_sorted.bam.bai")

        pysam.index(args.output)
