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
from pysam import reference

logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__ if __name__ != "__main__" else "rematch_reads_to_haplotypes")


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

def run_alignment(fa_seq, sequence, start_pos, sc_length, hap_cigar, aligner):
    """
    Perform the alignment between two sequences
    @param fa_seq: The reference sequence
    @param sequence: The query sequence
    @param start_pos: The start position of the read
    @param sc_length: The length of soft clipping
    @param aligner: The alignment object
    @return: The alignment score, CIGAR string, start position, and query start position
    """
    # Perform the alignment between two sequences
    if len(fa_seq) == 0 or len(sequence) == 0:
        return 0, 0, 0
    for alignment in aligner.align(fa_seq, sequence):
        # Print each alignment's score and the alignment itself
        start_pos_adjust, end_pos_adjust = adjust_start_end_positions(start_pos - sc_length, alignment.aligned, alignment.length, hap_cigar)

        return alignment.score, start_pos_adjust, end_pos_adjust
    return 0, 0, 0

def adjust_start_end_positions(start_pos, aligned, alignmnet_length, hap_cigar_tuples):
    """
        Convert the aligned segments in biopython format to a CIGAR farmat.
    """
    target_aligned, query_aligned = aligned
    adjusted_start_pos = 0  # This will store the adjusted start position based on initial target insertions

    if aligned.size == 0:
        return start_pos + adjusted_start_pos, start_pos + adjusted_start_pos


    t_gap = target_aligned[0][0]
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




def find_best_haplotype(read, local_aligner, reference):
    # in case cigar starts with soft clip, we need to adjust the haplotype start
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
    # extract the region of the assembly that the read overlaps
    with pysam.AlignmentFile(args.assembly, "rb") as assembly:
        for hap in assembly.fetch(read.reference_name, read_start, read_end):
            if not (hap.flag & 1024) and (
                    category == 1 or any(op in {1, 2, 4} and length > 20 for op, length in (hap.cigartuples or []))):
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
                score, start_pos_local, end_pos_local = run_alignment(hap_seq,
                                                                      read_seq,
                                                                      hap_start_position,
                                                                      0,
                                                                      hap.cigartuples,
                                                                      local_aligner)
                if (# overlap check
                    (read_start_position <= hap_end_position) and
                    (read_end_position >= hap_start_position) and
                    ((sc_size_start== 0 and sc_size_end == 0) or (sc_size_start > 0 and read.reference_start >= hap_start_position) or
                    (sc_size_end > 0 and read.reference_end <= hap_end_position))
                        and (score > best_score)):
                    best_score = score
                    best_hap = hap
                    best_start_point = start_pos_local - hap_start_position
                    best_end_point = end_pos_local - hap_start_position

        # check reference sequence as well
        # sift clip length from the start and end of the read
        sc_length_start = read.cigar[0][1] if read.cigar[0][0] == 4 else 0
        sc_length_end = read.cigar[-1][1] if read.cigar[-1][0] == 4 else 0
        del_length = sum([length for op, length in read.cigartuples if op == 2])

        ref_seq = reference[read.reference_name][
                 max(read.reference_start - sc_length_start, 0):
                 min(read.reference_end + sc_length_end + del_length, len(reference[read.reference_name]))].seq.upper()
        score, start_pos_local, end_pos_local = run_alignment(ref_seq,
                                                                read_seq,
                                                                read_start,
                                                                sc_size_start,
                                                                [],
                                                                local_aligner)
        if score > best_score:
            best_score = score
            best_hap = None
            best_start_point = start_pos_local - read_start
            best_end_point = end_pos_local - read_start


    # in case hap direction is reverse, we need to reverse the start and end points
    if best_hap is not None and best_hap.is_reverse:
        best_start_point, best_end_point = max(best_hap.query_length - best_end_point + 1, 0), best_hap.query_length - best_start_point + 1

    return best_hap, best_score, best_start_point, best_end_point




parser = argparse.ArgumentParser(description='Rematch reads to haplotypes')
parser.add_argument('--assembly', required=True, help='The input assembly file')
parser.add_argument('--reference', required=True, help='The reference genome FASTA file')
parser.add_argument('--output', required=True, help='The output assembly file with realigned supporting reads to haplotypes')
parser.add_argument("--n_jobs", help="n_jobs of parallel on contigs", type=int, default=-1)
parser.add_argument("--tumor_crams", help="The input tumor CRAM files", nargs='+', required=False)
parser.add_argument("--germline_crams", help="The input germline CRAM files", nargs='+', required=True)
parser.add_argument("--region", help="The region to process, in the format of chr<chr_num>:pos-pos", required=False)
args = parser.parse_args()

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

local_aligner = create_aligner('local', match, mismatch, gap_penalty, gap_extension_penalty, 0)
reference = pyfaidx.Fasta(args.reference, build_index=False)

haps_map = dict()
# align tumor and germline reads to haplotype
for cram_file, category in [(cram, 0) for cram in args.tumor_crams] + [(cram, 1) for cram in args.germline_crams]:
    with pysam.AlignmentFile(cram_file) as reads_cram:
        for read in reads_cram.fetch() if args.region is None else reads_cram.fetch(region = args.region):
            # Exclude PCR/optical duplicates, keep only tumor reads with long insertions, deletions, or soft-clips
            if not (read.flag & 1024) and (category == 1 or any(op in {1, 2, 4} and length > 20 for op, length in (read.cigartuples or []))):
                best_hap, best_score, start_point, end_point = find_best_haplotype(read, local_aligner, reference)
                if best_hap is not None:
                    if best_hap not in haps_map:
                        haps_map[best_hap] = []
                    # Append the SupportingRead object to the list
                    haps_map[best_hap].append(
                        SupportingRead(read.query_name, start_point, end_point, best_score, category))


# write the haplotypes to the output file
print("Writing the haplotypes to the output file")
with pysam.AlignmentFile(args.assembly) as assembly:
    with pysam.AlignmentFile(args.output + "_unsorted.bam", "wb", template=assembly) as output:
        # write each haplotype as a row in the output file
        # supporting reads are stored in the ef tag
        for hap, supporting_reads in haps_map.items():
            hap.set_tag("ef", " ".join([read.read_name for read in supporting_reads]))
            hap.set_tag("ez", " ".join([read.read_name for read in supporting_reads]))

            hap.set_tag("eq", array.array("f",[read.score for read in supporting_reads]))
            hap.set_tag("os", array.array("i",[read.start_overlap for read in supporting_reads]))
            hap.set_tag("oe", array.array("i",[read.end_overlap for read in supporting_reads]))
            hap.set_tag("ec", array.array("i",[read.category for read in supporting_reads]))
            hap.set_tag("et", array.array("b",[0 for read in supporting_reads])) # ?? Not sure about what is et
            output.write(hap)

# index output file
subprocess.check_call(f"samtools sort {args.output}_unsorted.bam -o {args.output}", shell=True)
subprocess.check_call(f"samtools index {args.output}", shell=True)

# remove unsoreted file
subprocess.check_call(f"rm {args.output}_unsorted.bam", shell=True)