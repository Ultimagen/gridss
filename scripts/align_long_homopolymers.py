# DESCRIPTION
#    This script realigns haplotypes in the areas that contain long homopolymer runs,
#    where UG data introduces false variation due to the limit on calling homopolymer length

import pysam
import pyfaidx
import tqdm.auto as tqdm
import argparse
import itertools
import re
from Bio import Align
from joblib import Parallel, delayed
import os
import logging
from Bio.Seq import Seq

logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__ if __name__ != "__main__" else "align_long_homopolymers")


def adjust_and_merge_cigar(cigar):
    """
    Adjusts for alternating deletion-insertion sequences, merges similar adjacent operations,
    and ensures preservation of terminal operations.
    @param cigar: The CIGAR string to adjust
    @return: The adjusted and merged CIGAR string
    """

    def calculate_adjustment(match):
        """
        in case of alternating deletion-insertion sequences, merge them
        """
        length1, op1, length2, op2 = match.groups()
        length1, length2 = int(length1), int(length2)

        if op1 == op2:
            return f"{length1 + length2}{op1}"
        else:
            overlap = min(length1, length2)
            leftover = abs(length1 - length2)
            adjusted_sequence = f'{overlap}M'
            if leftover > 0:
                leftover_op = op1 if length1 > length2 else op2
                adjusted_sequence += f'{leftover}{leftover_op}'
            return adjusted_sequence

    def merge_similar_operations(cigar):
        """
        Merges consecutive operations of the same type in a CIGAR string.
        """
        merged_cigar = []
        last_op = None
        last_length = 0

        for op, length in cigar_string_to_cigartuples(cigar, convert_int=False):
            length = int(length)
            if op == last_op:
                last_length += length
            else:
                if last_op is not None:
                    merged_cigar.append(f"{last_length}{last_op}")
                last_op = op
                last_length = length

        if last_op is not None:
            merged_cigar.append(f"{last_length}{last_op}")

        return ''.join(merged_cigar)

    adjusted_cigar = re.sub(r'(\d+)([DI])(\d+)([DI])', calculate_adjustment, cigar)
    merged_cigar = merge_similar_operations(adjusted_cigar)

    return merged_cigar


def cigar_may_change(cigartuples):
    """
    False only when adjust_and_merge_cigar would return the CIGAR unchanged:
    no two adjacent I/D operations, no two adjacent equal operations, no B operation (dropped by the parser)
    """
    if not cigartuples:
        return True
    prev_op = None
    for op, _ in cigartuples:
        if op == 9 or op == prev_op or (op in (1, 2) and prev_op in (1, 2)):
            return True
        prev_op = op
    return False


## ONLY FOR DEBUGGING
def calculate_sequence_length_by_cigar(cigar_string):
    """
    Calculate the length of the query sequence based on the CIGAR string
    """
    # Regular expression to find numbers followed by CIGAR operation characters
    pattern = re.compile(r'(\d+)([MIS=X])')
    seq_length = 0

    # Find all matches of the pattern in the CIGAR string
    matches = pattern.findall(cigar_string)

    # Sum up the lengths of operations that consume the query sequence
    for match in matches:
        count, operation = match
        if operation in "MIS=X":  # Check if the operation consumes the query sequence
            seq_length += int(count)

    return seq_length


cigar_op_codes = {
    'M': 0,  # BAM_CMATCH
    'I': 1,  # BAM_CINS
    'D': 2,  # BAM_CDEL
    'N': 3,  # BAM_CREF_SKIP
    'S': 4,  # BAM_CSOFT_CLIP
    'H': 5,  # BAM_CHARD_CLIP
    'P': 6,  # BAM_CPAD
    '=': 7,  # BAM_CEQUAL
    'X': 8  # BAM_CDIFF
}


# ONLY FOR DEBUGGING
def cigar_string_to_cigartuples(cigar_string, convert_int=True):
    """
    Convert a CIGAR string to a list of tuples of the form (operation, length)
    @param cigar_string: The CIGAR string to convert
    @param convert_int: Whether to convert the operation to an integer by cigar_op_codes (default: True)
    """
    # Regular expression to match each operation in the CIGAR string
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    cigartuples = []

    # Convert the operations to the format required by pysam
    if convert_int:
        for length, op in cigar_ops:
            cigartuples.append((cigar_op_codes[op], int(length)))
    else:
        # Do not convert the operation to an integer
        for length, op in cigar_ops:
            cigartuples.append((op, int(length)))

    return cigartuples


def remove_long_homopolymers(sequence, homopolymer_length=10):
    """
    Remove long homopolymers from the sequence and replace them with a homopolymer of length homopolymer_length
    @param sequence: The sequence to edit
    @param homopolymer_length: The length of homopolymers to search for
    @return: The edited sequence, the start and end positions of the deletions, and the total length of deletions
    """
    edited_parts = []
    start_end_del_tuples = []
    del_length = 0
    last_end = 0
    for run in re.finditer(r"(.)\1{%d,}" % homopolymer_length, sequence):
        start, end = run.span()
        edited_parts.append(sequence[last_end:start])
        edited_parts.append(sequence[start] * homopolymer_length)
        start_end_del_tuples.append((start + homopolymer_length, end))
        del_length += end - start - homopolymer_length
        last_end = end
    edited_parts.append(sequence[last_end:])
    return "".join(edited_parts), start_end_del_tuples, del_length


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


def alignment_score(fa_seq, sequence, aligner):
    """
    The score run_alignment would return, without building the alignment (traceback)
    """
    if len(fa_seq) == 0 or len(sequence) == 0:
        return 0
    return aligner.score(fa_seq, sequence)


def choose_global_or_local(fa_seq, sequence, ref_start, sc_length, global_aligner, local_aligner):
    """
    Choose between the global and the local alignment of the sequence, penalizing each soft clip of the local one
    The global alignment is built only when it is chosen
    @return: The score, CIGAR string, start position, reference start position, and whether global was chosen
    """
    score_local, cigar_local, start_pos_local, _, r_start_local = run_alignment(fa_seq, sequence, ref_start,
                                                                                sc_length, local_aligner)
    score_glob = alignment_score(fa_seq, sequence, global_aligner)
    penalize_times = cigar_local.count('S')
    if score_glob >= score_local + global_aligner.query_end_gap_score * penalize_times:
        _, cigar_glob, start_pos_glob, _, r_start_glob = run_alignment(fa_seq, sequence, ref_start, sc_length,
                                                                       global_aligner)
        return score_glob, cigar_glob, start_pos_glob, r_start_glob, True
    return score_local, cigar_local, start_pos_local, r_start_local, False


def align_and_choose(read, sequence, global_aligner, local_aligner, fa_seq, edited_sequence, ref_start, sc_length,
                     start_end_del_tuples):
    """
    Perform global/local alignment on the original sequence and the reverse complement, and choose the best alignment
    @param read: The read object
    @param sequence: The query sequence
    @param global_aligner: The global alignment object
    @param local_aligner: The local alignment object
    @param fa_seq: The reference sequence
    @param edited_sequence: The edited query sequence
    @param ref_start: The reference start position
    @param sc_length: The length of soft clipping
    @param start_end_del_tuples: The start and end positions of the deletions
    @return: The read object, the CIGAR string, the start position, the reference start position, the choice of alignment, and the start and end positions of the deletions
    the choice of alignment is 0 for the global alignment on the original sequence, 1 for the local alignment on the original sequence,
    2 for the global alignment on the reverse complement, and 3 for the local alignment on the reverse complement
    """
    score_orig, cigar, start_pos, r_start, orig_is_global = choose_global_or_local(fa_seq, edited_sequence, ref_start,
                                                                                    sc_length, global_aligner,
                                                                                    local_aligner)
    choice = 0 if orig_is_global else 1

    rev_comp_edited_sequence = str(Seq(edited_sequence).reverse_complement())
    # A local score is never below the global one (its end gaps cost <= 0), so it bounds the reverse score;
    # at or below it the original sequence wins (ties included) without the reverse alignments, the common case
    rev_score_bound = alignment_score(fa_seq, rev_comp_edited_sequence, local_aligner)
    if score_orig < rev_score_bound:
        score_rev, cigar_rev, start_pos_rev, r_start_rev, rev_is_global = choose_global_or_local(
            fa_seq, rev_comp_edited_sequence, ref_start, sc_length, global_aligner, local_aligner)
        if score_orig < score_rev:
            cigar = cigar_rev
            start_pos = start_pos_rev
            r_start = r_start_rev
            start_end_del_tuples = [(len(sequence) - end, len(sequence) - start) for start, end in start_end_del_tuples]
            start_end_del_tuples.reverse()
            choice = 2 if rev_is_global else 3

            read_qual = read.qual
            read.query_sequence = str(Seq(sequence).reverse_complement())
            read.qual = read_qual[::-1]

    return read, cigar, start_pos, r_start, choice, start_end_del_tuples


AMBIGUOUS_BASE = re.compile(r"[^ACGTacgt]")


def read_is_candidate_for_realign(read: pysam.AlignedSegment, min_required_softclip_length: int):
    '''Check if the read is a candidate for realignment (if it has minimal softclip of at least min_required_softclip_length)
    @param read: The read to check
    @param min_required_softclip_length: The minimum required soft clip length in the alignment to attempt realignment
    '''
    if read.cigartuples is None:
        return False
    if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] >= min_required_softclip_length:
        return True
    if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] >= min_required_softclip_length:
        return True
    return False

def realign_homopolymers(cram_path, output_path, reference_path, homopolymer_length=10, contig=None, min_required_softclip_length=10):
    """
    Find and realign homopolymers in reads in a CRAM file
    @param cram_path: The input CRAM file
    @param output_path: The output CRAM file
    @param reference_path: The reference genome FASTA file path
    @param homopolymer_length: The length of homopolymers to search for
    @param contig: The contig to process
    @param min_required_softclip_length: The minimum required soft clip length in the alignment to attempt realignment
    """
    # Open the CRAM file
    reference = pyfaidx.Fasta(reference_path, build_index=False)


    match = 1
    mismatch = -4
    gap_penalty = -6
    gap_extension_penalty = -1
    sc_penalty = -5

    global_aligner = create_aligner('global', match, mismatch, gap_penalty, gap_extension_penalty, sc_penalty)
    local_aligner = create_aligner('local', match, mismatch, gap_penalty, gap_extension_penalty, 0)

    # collect some statistics
    count = 0
    count_homopolymer = 0
    choices = [0, 0, 0, 0]

    with pysam.AlignmentFile(cram_path) as cram:
        with pysam.AlignmentFile(output_path, mode='wb', header=cram.header) as output:
            # Iterate over the reads in the CRAM file
            for read in cram.fetch(contig):
                count += 1
                sequence = read.query_sequence
                # Skip reads without a query sequence, or with ambiguous bases
                skip_read = sequence is None or not read_is_candidate_for_realign(read, min_required_softclip_length)
                skip_read = skip_read or AMBIGUOUS_BASE.search(sequence) is not None
                if not skip_read:
                    # Search for homopolymers in the sequence
                    edited_sequence, start_end_del_tuples, del_length = remove_long_homopolymers(read.query_sequence,
                                                                                                 homopolymer_length)
                    # Search for homopolymer in the reference
                    sc_length = 0
                    if read.cigartuples[0][0] == 4:
                        # the read is soft clipped
                        sc_length = read.cigartuples[0][1]
                    del_length = sum([length for op, length in read.cigartuples if op == 2])
                    fa_seq = reference[read.reference_name][
                             max(read.reference_start - sc_length, 0): read.reference_start + len(
                                 sequence) - sc_length + del_length].seq.upper()
                    edited_fa_seq, ref_start_end_del_tuples, ref_del_length = remove_long_homopolymers(fa_seq,
                                                                                                       homopolymer_length)

                    if ref_del_length > del_length:
                        # In case the reference has more deletions than the read, we need to adjust the read
                        fa_seq = reference[read.reference_name][
                                 max(read.reference_start - sc_length, 0): read.reference_start + len(
                                     sequence) - sc_length + ref_del_length - del_length].seq.upper()
                        edited_fa_seq, ref_start_end_del_tuples, ref_del_length = remove_long_homopolymers(fa_seq,
                                                                                                           homopolymer_length)

                if not skip_read and (len(start_end_del_tuples) > 0 or len(ref_start_end_del_tuples) > 0):
                    # we found long homopolymer and want to run alignment on that with homopolymere of the length of homopolymer_length
                    count_homopolymer += 1
                    read, cigar, start_pos, r_start, choice, start_end_del_tuples = align_and_choose(read, sequence,
                                                                                                     global_aligner,
                                                                                                     local_aligner,
                                                                                                     edited_fa_seq,
                                                                                                     edited_sequence,
                                                                                                     read.reference_start,
                                                                                                     sc_length,
                                                                                                     start_end_del_tuples)
                    choices[choice] += 1

                    if choice == 0 | choice == 1:
                        read.is_reverse = False
                    else:
                        read.is_reverse = True

                    # Add insertions and deletions into the CIGAR string
                    updated_cigar = cigar
                    # Insertions by the query
                    for start_del, end_del in start_end_del_tuples:
                        updated_cigar, _ = insert_operation_into_cigar(updated_cigar, start_del, end_del - start_del, 0,
                                                                       'I')
                    # Deletions by the reference
                    for start_del, end_del in ref_start_end_del_tuples:
                        updated_cigar, start_pos = insert_operation_into_cigar(updated_cigar, start_del - r_start,
                                                                               end_del - start_del, start_pos, 'D')

                    adjusted_cigar = adjust_and_merge_cigar(updated_cigar)
                    read.cigar = cigar_string_to_cigartuples(adjusted_cigar)
                    # edge case when the start position is negative
                    # Can happen when the read is soft-clipped at the beginning of the chromosome
                    if read.reference_start - sc_length < 0:
                        read.reference_start = start_pos - (read.reference_start - sc_length)
                    else:
                        read.reference_start = start_pos

                    if not (len(sequence) == calculate_sequence_length_by_cigar(
                            adjusted_cigar) or calculate_sequence_length_by_cigar(adjusted_cigar) == 0):
                        raise AssertionError(
                            f"failed for read: {read.query_name}\ncigar: {adjusted_cigar}\ncigar_len:{calculate_sequence_length_by_cigar(adjusted_cigar)}\nseq_len:{len(sequence)}\nsequence: {sequence}\nedited_sequence: {edited_sequence}")

                    output.write(read)
                else:
                    if cigar_may_change(read.cigartuples):
                        adjusted_cigar = adjust_and_merge_cigar(read.cigarstring)
                        if adjusted_cigar != read.cigarstring:
                            read.cigar = cigar_string_to_cigartuples(adjusted_cigar)
                    output.write(read)

    reference.close()
    return count, count_homopolymer, choices



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
        return 0, "", 0, 0, 0
    for alignment in aligner.align(fa_seq, sequence):
        # Print each alignment's score and the alignment itself
        start_pos_adjust, cigar, q_start, r_start = convert_alignment_to_cigar(alignment.aligned, len(sequence))
        start_pos = start_pos - sc_length + start_pos_adjust
        return alignment.score, cigar, start_pos, q_start, r_start
    return 0, "", 0, 0, 0


def insert_operation_into_cigar(cigar, position, op_size, start_pos, op_type):
    """
    Insert an operation into a CIGAR string at a specified position.
    @param cigar: The CIGAR string to modify
    @param position: The position at which to insert the operation
    @param op_size: The size of the operation to insert
    @param start_pos: The start position of the read
    @param op_type: The type of operation to insert ('D' for deletion or 'I' for insertion)
    @return: The updated CIGAR string and the updated start position
    """
    assert op_type in ['D', 'I'], "op_type must be 'D' for deletion or 'I' for insertion"

    # Split the CIGAR string into its components
    parsed_cigar = cigar_string_to_cigartuples(cigar, convert_int=False)
    new_cigar_ops = []
    accumulated_length = 0
    operation_inserted = False
    i = 0
    if position < 0 and op_type == 'D':
        start_pos += op_size
        return cigar, start_pos
    while i < len(parsed_cigar):
        op, count = parsed_cigar[i]
        skip_operation = False
        if op_type == 'D':  # Increment accumulated_length except for 'D' operation in when the position is by the query
            if op == 'I' or op == 'S':
                skip_operation = True
        elif op_type == 'I':  # Increment accumulated_length except for 'I' operation in when the position is by the reference
            if op == 'D':
                skip_operation = True

        if op in 'MIDS' and not operation_inserted:
            if (not skip_operation) and (accumulated_length <= position <= accumulated_length + count):
                # Split the operation at the insertion/deletion position
                operation_inserted = True
                if op == op_type:
                    # merge the operations
                    new_cigar_ops.append((op_size + count, op_type))
                elif op == 'S':
                    # merge the operations
                    new_cigar_ops.append((op_size + count, 'S'))
                else:
                    # Insert the operation
                    before_pos = position - accumulated_length
                    if before_pos > 0:
                        new_cigar_ops.append((before_pos, op))

                    after_pos = count - before_pos
                    if after_pos == 0:
                        # check if we can merge with the next operation
                        if i + 1 < len(parsed_cigar) and parsed_cigar[i + 1][1] == op_type:
                            # merge the operations
                            new_cigar_ops.append((op_size + parsed_cigar[i + 1][0], op_type))
                            i = i + 1
                        else:
                            new_cigar_ops.append((op_size, op_type))
                    elif after_pos > 0:
                        new_cigar_ops.append((op_size, op_type))
                        new_cigar_ops.append((after_pos, op))

            else:
                new_cigar_ops.append((count, op))

            if not skip_operation:
                accumulated_length += count

        else:
            new_cigar_ops.append((count, op))
        i = i + 1

    # Convert back to CIGAR string format
    updated_cigar = ''.join(f"{count}{op}" for count, op in new_cigar_ops)

    return updated_cigar, start_pos


def realign_and_sort(input_path, part_path, reference_path, homopolymer_length, contig,
                     min_required_softclip_length):
    """
    realign_homopolymers on one contig, then coordinate sort its output into part_path
    """
    unsorted_path = f"{part_path}.unsorted.bam"
    result = realign_homopolymers(input_path, unsorted_path, reference_path, homopolymer_length, contig,
                                  min_required_softclip_length)
    pysam.sort("--no-PG", "-o", part_path, unsorted_path)
    os.remove(unsorted_path)
    return result


def copy_contigs(input_path, part_path, contigs):
    """
    Copy the reads of short contigs unchanged; returns the same statistics as realign_homopolymers
    """
    count = 0
    with pysam.AlignmentFile(input_path) as bam:
        with pysam.AlignmentFile(part_path, mode='wb', header=bam.header) as output:
            for contig in contigs:
                for read in bam.fetch(contig):
                    output.write(read)
                    count += 1
    return count, 0, [0, 0, 0, 0]


def main():
    parser = argparse.ArgumentParser(description='Find and realign homopolymers in reads.')
    parser.add_argument('--input', required=True, help='The input CRAM file')
    parser.add_argument('--output', required=True, help='The output CRAM file')
    parser.add_argument('--reference', required=True, help='The reference genome FASTA file')
    parser.add_argument('--homopolymer_length', type=int, default=10, help='The length of homopolymers to search for')
    parser.add_argument('--min_required_softclip_length', type=int, default=10, help='The minimum required soft clip length in the alignment to attempt realignment')
    parser.add_argument("--n_jobs", help="n_jobs of parallel on contigs", type=int, default=-1)
    args = parser.parse_args()

    MIN_CONTIG_LENGTH = 100000

    with pysam.AlignmentFile(args.input, "rc") as cram_file:
        contig_lengths = dict(zip(cram_file.references, cram_file.lengths))
        contig_reads = {stat.contig: stat.total for stat in cram_file.get_index_statistics()}

    # One part file per large contig and per run of short contigs, in header order; concatenated in that order
    parts = []
    for is_large, group in itertools.groupby(contig_lengths, key=lambda c: contig_lengths[c] > MIN_CONTIG_LENGTH):
        group = list(group)
        if is_large:
            parts.extend((f"{args.output}{contig}.bam", [contig], True) for contig in group)
        else:
            parts.append((f"{args.output}{group[0]}.short_contigs.bam", group, False))

    # Largest first, so that no large contig starts last and leaves the other workers idle
    by_size = sorted(parts, key=lambda part: -sum(contig_reads.get(contig, 0) for contig in part[1]))
    results = Parallel(n_jobs=args.n_jobs, backend="multiprocessing", max_nbytes=None)(
        delayed(realign_and_sort)(args.input, part_path, args.reference, args.homopolymer_length, contigs[0],
                                  args.min_required_softclip_length)
        if is_large else delayed(copy_contigs)(args.input, part_path, contigs)
        for part_path, contigs, is_large in tqdm.tqdm(by_size)
    )
    total_count = 0
    total_count_homopolymer = 0
    total_choices = [0, 0, 0, 0]
    for count, count_homopolymer, choices in results:
        total_count += count
        total_count_homopolymer += count_homopolymer
        total_choices = [total_choices[i] + choices[i] for i in range(4)]

    part_paths = [part_path for part_path, _, _ in parts]
    pysam.cat("--no-PG", "-o", args.output, *part_paths)
    for part_path in part_paths:
        os.remove(part_path)
    pysam.index(args.output)

    logger.info(f"Total reads: {total_count}")
    logger.info(f"Reads with homopolymers: {total_count_homopolymer}")
    logger.info(f"Reads where original sequence is better: {total_choices[0]}")
    logger.info(f"Reads where original sequence is better (local): {total_choices[1]}")
    logger.info(f"Reads where reverse complement sequence is better: {total_choices[2]}")
    logger.info(f"Reads where reverse complement sequence is better (local): {total_choices[3]}")


if __name__ == "__main__":
    main()
