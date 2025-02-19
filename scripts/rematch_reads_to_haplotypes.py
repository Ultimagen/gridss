# DESCRIPTION
#    This script realigns haplotypes in the areas that contain long homopolymer runs,
#    where UG data introduces false variation due to the limit on calling homopolymer length
import numpy as np
import pysam
import argparse
import logging
import subprocess
from Bio import Align
import array

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

def run_alignment(fa_seq, sequence, start_pos, sc_length, aligner):
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
        return 0, "", 0, 0, 0, 0
    for alignment in aligner.align(fa_seq, sequence):
        # Print each alignment's score and the alignment itself
        start_pos_adjust, cigar, q_start, r_start = convert_alignment_to_cigar(alignment.aligned, len(sequence))
        start_pos = start_pos - sc_length + start_pos_adjust
        return alignment.score, cigar, start_pos, q_start, r_start, alignment.length
    return 0, "", 0, 0, 0, 0

def convert_alignment_to_cigar(aligned, seq2_len):
    """
        Convert the aligned segments in biopython format to a CIGAR farmat.
    """
    target_aligned, query_aligned = aligned
    cigar = []
    adjusted_start_pos = 0  # This will store the adjusted start position based on initial target insertions

    if aligned.size == 0:
        return adjusted_start_pos, '0M', 0, 0

    # Helper function to add operation to CIGAR
    def add_op(op, length):
        if length > 0:
            cigar.append(f"{length}{op}")

    last_target_end, last_query_end = 0, 0
    first_time = True
    for (t_start, t_end), (q_start, q_end) in zip(target_aligned, query_aligned):
        # Handle gaps before this aligned segment
        t_gap = t_start - last_target_end
        q_gap = q_start - last_query_end

        # Determine initial soft clipping for the query and adjust start_pos for target insertions
        if first_time:
            first_time = False
            if t_gap > 0:
                adjusted_start_pos = t_gap
            if q_gap > 0:
                add_op('S', q_gap)

        # Adjust for gaps
        # Compare t_gap and q_gap to determine order
        elif t_gap > 0 and q_gap > 0:
            # Both gaps exist, determine order
            if t_start < q_start:
                # Deletion occurs first
                add_op('D', t_gap)
                add_op('I', q_gap)
            else:
                # Insertion occurs first
                add_op('I', q_gap)
                add_op('D', t_gap)
        else:
            # Only one type of gap exists
            if t_gap > 0: add_op('D', t_gap)
            if q_gap > 0: add_op('I', q_gap)

        # Add matched segment
        segment_length = min(t_end - t_start, q_end - q_start)
        add_op('M', segment_length)

        # Update last processed positions
        last_target_end, last_query_end = t_end, q_end

    # Handle gaps after the last aligned segment - add soft-clipping for the remaining query
    if last_query_end < seq2_len:
        add_op('S', seq2_len - last_query_end)  # Treat remaining query as soft clipped

    return adjusted_start_pos, ''.join(cigar), query_aligned[0][0], target_aligned[0][0]  # , query_aligned[-1][1]

def find_best_haplotype(read, local_aligner):
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

    # extract the region of the assembly that the read overlaps
    subprocess.check_call(
        f"samtools view -b {args.assembly} {read.reference_name}:{read_start}-{read_end} -o {args.output}_temp_hap.bam",
        shell=True)
    subprocess.check_call(f"samtools index {args.output}_temp_hap.bam", shell=True)
    # open the extracted region and find the best haplotype
    with (pysam.AlignmentFile(args.output + "_temp_hap.bam") as hap_cram):
        # find the best alignment
        best_score = -np.inf
        best_hap = None
        best_start_point = None
        best_end_point = None
        read_seq = read.query_sequence
        for hap in hap_cram.fetch():


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
            score, cigar, start_pos_local, q_start_local, r_start_local, align_length = run_alignment(hap_seq,
                                                                                              read_seq,
                                                                                              hap_start_position,
                                                                                              0,
                                                                                              local_aligner)

            print(
                f"Processing haplotype {hap.query_name} Alignment score: {score} start_pos_local: {start_pos_local} end position: {start_pos_local + align_length}")

            if (# overlap check
                (read_start_position <= hap_end_position) and
                (read_end_position >= hap_start_position) and
                ((sc_size_start== 0 and sc_size_end == 0) or (sc_size_start > 0 and read.reference_start >= hap_start_position) or
                (sc_size_end > 0 and read.reference_end <= hap_end_position))
                    and (score > best_score)):
                best_score = score
                best_hap = hap
                best_start_point = start_pos_local - hap_start_position
                best_end_point = start_pos_local + align_length - hap_start_position
        if best_hap is not None:
            print(f"Best haplotype for read {read.query_name} is {best_hap.query_name} with score {best_score} start point {best_start_point} end point {best_end_point}")
        else:
            print(f"No haplotype found for read {read.query_name}")
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



if(args.region is not None):

    logger.info(f"Processing region {args.region}")
    # assembly file is in CRAM format, extract the region
    subprocess.check_call(
        f"samtools view -b {args.assembly} {args.region} -o {args.output}_temp_{args.region}.bam",
        shell=True)
    subprocess.check_call(f"samtools index {args.output}_temp_{args.region}.bam", shell=True)
    args.assembly = args.output + f"_temp_{args.region}.bam"
    # tumor crams, extract the region
    for i in range(len(args.tumor_crams)):
        subprocess.check_call(
            f"samtools view -b {args.tumor_crams[i]} {args.region} -o {args.output}_temp_tumor_{args.region}_{i}.bam",
            shell=True)
        subprocess.check_call(f"samtools index {args.output}_temp_tumor_{args.region}_{i}.bam", shell=True)
        args.tumor_crams[i] = args.output + f"_temp_tumor_{args.region}_{i}.bam"

    # germline crams, extract the region
    for i in range(len(args.germline_crams)):
        subprocess.check_call(
            f"samtools view -b {args.germline_crams[i]} {args.region} -o {args.output}_temp_germline_{args.region}_{i}.bam",
            shell=True)
        subprocess.check_call(f"samtools index {args.output}_temp_germline_{args.region}_{i}.bam", shell=True)
        args.germline_crams[i] = args.output + f"_temp_germline_{args.region}_{i}.bam"

    tumor_crams = [f"{args.output}_temp_germline_{args.region}_{i}.bam" for i in range(len(args.tumor_crams))]
    germline_crams = [f"{args.output}_temp_germline_{args.region}_{i}.bam" for i in range(len(args.germline_crams))]
    merged_tumor_cram = f"{args.output}_temp_merged_tumor_{args.region}.bam"
    merged_germline_cram = f"{args.output}_temp_merged_germline_{args.region}.bam"
else:
    tumor_crams = args.tumor_crams
    germline_crams = args.germline_crams
    merged_tumor_cram = f"{args.output}_temp_merged_tumor.bam"
    merged_germline_cram = f"{args.output}_temp_merged_germline.bam"


if len(args.tumor_crams) > 0:
    subprocess.check_call(f"samtools merge -f {merged_tumor_cram} {' '.join([cram for cram in tumor_crams])}", shell=True)
    subprocess.check_call(f"samtools index {merged_tumor_cram}", shell=True)

subprocess.check_call(f"samtools merge -f {merged_germline_cram} {' '.join([cram for cram in germline_crams])}", shell=True)
subprocess.check_call(f"samtools index {merged_germline_cram}", shell=True)

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



haps_map = dict()

# align tumor reads to haplotype
if len(args.tumor_crams) > 0:
    with pysam.AlignmentFile(args.output+"_temp_merged_tumor.cram") as reads_cram:
            for read in reads_cram.fetch():
                print(f"Processing read {read.query_name} {read.reference_name}:{read.reference_start}-{read.reference_end}")
                best_hap, best_score, start_point, end_point = find_best_haplotype(read, local_aligner)
                if best_hap is not None:
                    if best_hap not in haps_map:
                        haps_map[best_hap] = []
                    # Append the SupportingRead object to the list
                    haps_map[best_hap].append(
                        SupportingRead(read.query_name, start_point, end_point, best_score, 0))

# align germline reads to haplotype
with pysam.AlignmentFile(args.output+"_temp_merged_germline.cram") as reads_cram:
        for read in reads_cram.fetch():
            print(f"Processing read {read.query_name} {read.reference_name}:{read.reference_start}-{read.reference_end}")
            best_hap, best_score, start_point, end_point = find_best_haplotype(read, local_aligner)
            if best_hap is not None:
                if best_hap not in haps_map:
                    haps_map[best_hap] = []
                # Append the SupportingRead object to the list
                haps_map[best_hap].append(
                    SupportingRead(read.query_name, start_point, end_point, best_score, 1))

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