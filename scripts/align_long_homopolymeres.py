import pysam
import pyfaidx
import argparse
import re
from Bio import Align
from joblib import Parallel, delayed
import os




def merge_similar_operations(cigar):
    """
    Merges consecutive operations of the same type.
    """
    merged_cigar = []
    last_op = None
    last_length = 0

    # Extract operations and their lengths
    for length, op in re.findall(r'(\d+)([MID])', cigar):
        length = int(length)
        if op == last_op:
            last_length += length
        else:
            if last_op is not None:
                merged_cigar.append(f"{last_length}{last_op}")
            last_length = length
            last_op = op

    # Don't forget to add the last operation
    if last_op is not None:
        merged_cigar.append(f"{last_length}{last_op}")

    return ''.join(merged_cigar)


def merge_similar_operations(cigar):
    """
    Merges consecutive operations of the same type in a CIGAR string.
    """
    merged_cigar = []
    last_op = None
    last_length = 0

    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
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

def adjust_and_merge_cigar(cigar):
    """
    Adjusts for alternating deletion-insertion sequences, merges similar adjacent operations,
    and ensures preservation of terminal operations.
    """
    def calculate_adjustment(match):
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

    adjusted_cigar = re.sub(r'(\d+)([DI])(\d+)([DI])', calculate_adjustment, cigar)
    merged_cigar = merge_similar_operations(adjusted_cigar)

    return merged_cigar


## ONLY FOR DEBUGGING
def calculate_sequence_length_by_cigar(cigar_string):
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
    'X': 8   # BAM_CDIFF
}
# ONLY FOR DEBUGGING
def cigar_string_to_cigartuples(cigar_string):
    # Regular expression to match each operation in the CIGAR string
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    cigartuples = []

    for length, op in cigar_ops:
        cigartuples.append((cigar_op_codes[op], int(length)))

    return cigartuples



def remove_long_homopolymers(sequence, homopolymer_length=10):
    i = 0
    del_length = 0
    edited_sequence = ""
    start_end_del_tuples = []
    while i < len(sequence):
        count = 1
        while i + count < len(sequence) and sequence[i] == sequence[i + count]:
            count += 1

        if count > homopolymer_length:
            edited_sequence += sequence[i] * homopolymer_length
            i += count
            start_end_del_tuples.append((i - (count - homopolymer_length) , i))
            del_length += count - homopolymer_length
        elif count > 1:
            edited_sequence += sequence[i] * count
            i += count
        else:
            edited_sequence += sequence[i]
            i += 1
    return edited_sequence, start_end_del_tuples, del_length


def find_homopolymers(cram_path, output_path, reference_path, homopolymer_length=10, contig = None):
    # Open the CRAM file
    reference = pyfaidx.Fasta(reference_path)

    match = 1
    mismatch = -4
    gap_penalty = -6
    gap_extension_penalty = -1
    sc_penalty = -5

    # Create a new PairwiseAligner object
    global_aligner = Align.PairwiseAligner()
    global_aligner.mode = 'global'
    global_aligner.match_score = match
    global_aligner.mismatch_score = mismatch
    # Set gap scoring
    global_aligner.open_gap_score = gap_penalty
    global_aligner.extend_gap_score = gap_extension_penalty

    # Set end gap scores to 0.0 to not penalize them
    global_aligner.target_end_gap_score = 0
    global_aligner.query_end_gap_score = sc_penalty

    local_aligner = Align.PairwiseAligner()
    local_aligner.mode = 'local'
    local_aligner.match_score = match
    local_aligner.mismatch_score = mismatch
    # Set gap scoring
    local_aligner.open_gap_score = gap_penalty
    local_aligner.extend_gap_score = gap_extension_penalty
    local_aligner.target_end_gap_score = 0
    local_aligner.query_end_gap_score = 0


    # collect some statistics
    count = 0
    count_homopolymere  = 0
    count_orig_local = 0
    count_orig_glob = 0
    count_rev_local = 0
    count_rev_glob = 0


    with pysam.AlignmentFile(cram_path) as cram:
        with pysam.AlignmentFile(output_path, mode='wb', header=cram.header) as output:
            # Iterate over the reads in the CRAM file

            for read in cram.fetch(contig):
                count += 1
                sequence = read.query_sequence

                if sequence is None or 'Y' in sequence:
                    continue  # Skip reads without a query sequence, or with ambiguous bases

                # Search for homopolymers in the sequence
                edited_sequence, start_end_del_tuples, del_length = remove_long_homopolymers(read.query_sequence, homopolymer_length)

                #Search for homopolymer in the reference
                sc_length = 0
                if read.cigartuples[0][0] == 4:
                    # the read is soft clipped
                    sc_length = read.cigartuples[0][1]
                fa_seq = reference[read.reference_name][read.reference_start - sc_length - 0 : read.reference_start + len(sequence) - sc_length].seq.upper()
                edited_fa_seq, ref_start_end_del_tuples, ref_del_length = remove_long_homopolymers(fa_seq, homopolymer_length)

                if ref_del_length > del_length:

                    fa_seq = reference[read.reference_name][
                             read.reference_start - sc_length: read.reference_start + len(
                                 sequence) - sc_length + ref_del_length - del_length].seq.upper()
                    edited_fa_seq, ref_start_end_del_tuples, ref_del_length = remove_long_homopolymers(fa_seq,homopolymer_length)


                if len(start_end_del_tuples)>0 or len(ref_start_end_del_tuples)>0:
                    # we found long homopolymer and want to run alignment on that with homopolymere of the length of homopolymer_length
                    count_homopolymere += 1
                    # print(
                    #     f"Read {read.query_name}")
                    # print(start_end_del_tuples)
                    #
                    # print(f"Original sequence:  {sequence}")
                    # print(f"Edited sequence:    {edited_sequence}")
                    #
                    # print(f"Original reference: {fa_seq}")
                    # print(f"Edited reference:   {edited_fa_seq}")
                    #
                    # print("Read cigar: ", read.cigarstring)
                    #
                    # print("Running alignment on original sequence")
                    # print("Running global alignment")
                    #
                    # # run alignment on original sequence
                    # print("Running alignment on original sequence")
                    # print("Running global alignment")
                    score_orig_glob, cigar_orig_glob, start_pos_glob, q_start_glob, r_start_glob = run_alignment_biopyhon(edited_fa_seq, edited_sequence, read.reference_start, sc_length, 0, global_aligner)
                    # print("Running local alignment")
                    score_orig_local, cigar_orig_local, start_pos_local, q_start_local, r_start_local = run_alignment_biopyhon(edited_fa_seq, edited_sequence, read.reference_start, sc_length, 0, local_aligner)
                    penalize_times = cigar_orig_local.count('S')
                    if score_orig_glob >= score_orig_local + sc_penalty * penalize_times:
                        score_orig = score_orig_glob
                        cigar_orig = cigar_orig_glob
                        start_pos = start_pos_glob
                        r_start = r_start_glob
                        q_start = q_start_glob
                        chosen_orig = 1
                    else:
                        score_orig = score_orig_local
                        cigar_orig = cigar_orig_local
                        start_pos = start_pos_local
                        r_start = r_start_local
                        q_start = q_start_local
                        chosen_orig = 2

                    # run alignment on reverse complement

                    rev_comp_edited_sequence = reverse_complement(edited_sequence)
                    # print("Running alignment on reverse complement")
                    # print("Running global rev alignment")
                    score_rev_glob, cigar_rev_glob, start_pos_rev_glob, q_start_rev_glob, r_start_rev_glob = run_alignment_biopyhon(edited_fa_seq, rev_comp_edited_sequence, read.reference_start, sc_length, 0, global_aligner)
                    # print("Running local rev alignment")
                    score_rev_local, cigar_rev_local, start_pos_rev_local,  q_start_rev_local, r_start_rev_local = run_alignment_biopyhon(edited_fa_seq, rev_comp_edited_sequence,read.reference_start, sc_length, 0, local_aligner)

                    penalize_times = cigar_rev_local.count('S')
                    if score_rev_glob >= score_rev_local + sc_penalty * penalize_times:
                        score_rev = score_rev_glob
                        cigar_rev = cigar_rev_glob
                        start_pos_rev = start_pos_rev_glob
                        q_start_rev = q_start_rev_glob
                        r_start_rev = r_start_rev_glob
                        chosen_rev = 1
                    else:
                        score_rev = score_rev_local
                        cigar_rev = cigar_rev_local
                        start_pos_rev = start_pos_rev_local
                        q_start_rev = q_start_rev_local
                        r_start_rev = r_start_rev_local
                        chosen_rev = 2

                    if score_orig >= score_rev:

                        cigar = cigar_orig
                        if chosen_orig == 1:
                            count_orig_glob += 1
                        else:
                            count_orig_local += 1
                    else:
                        print(
                             f"rev Read {read.query_name}")
                        score = score_rev
                        cigar = cigar_rev
                        edited_seq = rev_comp_edited_sequence
                        start_pos = start_pos_rev
                        q_start = q_start_rev
                        r_start = r_start_rev
                        start_end_del_tuples = [(len(sequence) - end, len(sequence) - start) for start, end in start_end_del_tuples]
                        start_end_del_tuples.reverse()
                        #ref_start_end_del_tuples = [(len(fa_seq) - end, len(fa_seq) - start) for start, end in ref_start_end_del_tuples]
                        #ref_start_end_del_tuples.reverse()
                        if chosen_rev == 1:
                            count_rev_glob += 1
                            # print("rev_glob")
                        else:
                            count_rev_local += 1
                            # print("rev_local")
                        read.query_sequence  = reverse_complement(sequence)


                    # Insert deletions into the CIGAR string
                    updated_cigar = cigar

                    # print(f"Final cigar before change:        {cigar}")
                    # print("Final cigar before change (length): ", calculate_sequence_length_by_cigar(cigar))
                    # print(f"seq start_end_del_tuples: {start_end_del_tuples}")
                    # print(f"                             {updated_cigar}")

                    for start_del, end_del in start_end_del_tuples:
                        updated_cigar, _ = insert_operation_into_cigar(updated_cigar, start_del, end_del - start_del, 0, 'I')

                    updated_cigar_D = updated_cigar
                    # print(f"ref start_end_del_tuples: {ref_start_end_del_tuples}")
                    # print(f"                             {updated_cigar}")

                    for start_del, end_del in ref_start_end_del_tuples:
                        updated_cigar, start_pos = insert_operation_into_cigar(updated_cigar, start_del - r_start, end_del - start_del, start_pos, 'D')

                    # print(f"updated_cigar_D                   {updated_cigar_D}")
                    # print("Final cigar before change_D (length): ", calculate_sequence_length_by_cigar(updated_cigar_D))
                    # print(f"Final cigar after change:         {updated_cigar}")
                    # print("Final cigar after change (length): ", calculate_sequence_length_by_cigar(updated_cigar))
                    # print("cigartuples: ", cigar_string_to_cigartuples(updated_cigar))
                    adjusted_cigar = adjust_and_merge_cigar(updated_cigar)
                    read.cigar = cigar_string_to_cigartuples(adjusted_cigar)
                    read.reference_start = start_pos


                    assert (len(edited_sequence) == calculate_sequence_length_by_cigar(cigar_orig_local) or calculate_sequence_length_by_cigar(cigar_orig_local) == 0)
                    assert(len(edited_sequence) == calculate_sequence_length_by_cigar(cigar_orig_glob) or calculate_sequence_length_by_cigar(cigar_orig_glob) == 0)
                    assert (len(edited_sequence) == calculate_sequence_length_by_cigar(cigar_rev_local) or calculate_sequence_length_by_cigar(cigar_rev_local) == 0)
                    assert(len(edited_sequence) == calculate_sequence_length_by_cigar(cigar_rev_glob) or calculate_sequence_length_by_cigar(cigar_rev_glob) == 0)
                    if not (len(sequence) == calculate_sequence_length_by_cigar(updated_cigar) or calculate_sequence_length_by_cigar(updated_cigar) == 0):
                        raise AssertionError(f"failed for read: {read.query_name}\ncigar: {updated_cigar}\ncigar_len:{calculate_sequence_length_by_cigar(updated_cigar)}\nseq_len:{len(sequence)}\nsequence: {sequence}\nedited_sequence: {edited_sequence}")
                    if not (len(sequence) == calculate_sequence_length_by_cigar(adjusted_cigar) or calculate_sequence_length_by_cigar(adjusted_cigar) == 0):
                        raise AssertionError(f"failed for read: {read.query_name}\ncigar: {adjusted_cigar}\ncigar_len:{calculate_sequence_length_by_cigar(adjusted_cigar)}\nseq_len:{len(sequence)}\nsequence: {sequence}\nedited_sequence: {edited_sequence}")

                    print(
                        f"Read {read.query_name} {updated_cigar}")
                    output.write(read)
                else:
                    output.write(read)
                # print("###########")

    print(f"Total reads: {count}")
    print(f"Total homopolymer reads: {count_homopolymere}")
    print(f"Total reads where original sequence is better: {count_orig_glob}")
    print(f"Total reads where original sequence is better (local): {count_orig_local}")
    print(f"Total reads where reverse complement sequence is better: {count_rev_glob}")
    print(f"Total reads where reverse complement sequence is better (local): {count_rev_local}")



def convert_aligned_to_cigar(aligned, seq2_len):
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

    return adjusted_start_pos, ''.join(cigar), query_aligned[0][0], target_aligned[0][0]#, query_aligned[-1][1]

def run_alignment_biopyhon(fa_seq, sequence, start_pos, sc_length, del_length, aligner):
    # Perform the alignment between two sequences
    if len(fa_seq) == 0 or len(sequence) == 0:
        return 0, "", 0, 0, 0
    for alignment in aligner.align(fa_seq, sequence ):

        # Print each alignment's score and the alignment itself
        start_pos_adjust, cigar, q_start, r_start = convert_aligned_to_cigar(alignment.aligned,  len(sequence))

        start_pos = start_pos - sc_length - del_length + start_pos_adjust

        # print(f"Cigar = {cigar}")
        # print(alignment)
        return alignment.score, cigar, start_pos, q_start, r_start
    return 0, "", 0, 0, 0


def reverse_complement(seq):
    """Return the reverse complement of the given DNA sequence."""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    # Generate the complement sequence
    complement_seq = ''.join(complement_dict[nuc] for nuc in seq)
    # Reverse the complement sequence
    reverse_comp_seq = complement_seq[::-1]
    return reverse_comp_seq



def insert_operation_into_cigar(cigar, position, op_size, start_pos, op_type):
    assert op_type in ['D', 'I'], "op_type must be 'D' for deletion or 'I' for insertion"

    # Split the CIGAR string into its components
    cigar_components = re.findall(r'\d+[MIDNSHP=X]', cigar)

    # Convert components into a list of tuples (count, operation)
    parsed_cigar = [(int(comp[:-1]), comp[-1]) for comp in cigar_components]

    new_cigar_ops = []
    accumulated_length = 0
    operation_inserted = False
    i=0
    if position < 0 and op_type == 'D':
        start_pos += op_size
        return cigar, start_pos
    while i<len(parsed_cigar):
        count, op = parsed_cigar[i]
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
                    new_cigar_ops.append((op_size+count, op_type))
                elif op == 'S':
                    # merge the operations
                    new_cigar_ops.append((op_size+count, 'S'))
                else:
                    # Insert the operation
                    before_pos = position - accumulated_length
                    if before_pos > 0:
                        new_cigar_ops.append((before_pos, op))

                    after_pos = count - before_pos
                    if after_pos == 0:
                        # check if we can merge with the next operation
                        if i+1 < len(parsed_cigar) and parsed_cigar[i+1][1] == op_type:
                            # merge the operations
                            new_cigar_ops.append((op_size+parsed_cigar[i+1][0], op_type))
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
        i=i+1

    # Merge operations if needed
    # merged_ops = merge_operations(new_cigar_ops)

    # Convert back to CIGAR string format
    updated_cigar = ''.join(f"{count}{op}" for count, op in new_cigar_ops)

    return updated_cigar, start_pos

# # usage
# cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9_69449415_69450719.bam"
# output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9_69449415_69450719.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)
#

# cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9.bam"
# output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)


# cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9_6790002_6790202.bam"
# output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9_6790002_6790202.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)




# cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9_275409_339530.bam"
# output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9_275409_339530.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)

# cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9_1_275409.bam"
# output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9_1_275409.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)


# cram_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_chr1_94734178_94934178.bam"
# output_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_long_homopolymeres_aligned_unsorted_chr1_94734178_94934178.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)

# cram_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_chr1.bam"
# output_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_long_homopolymeres_aligned_unsorted_chr1.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)

# cram_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_chr1_248745192_248765192.bam"
# output_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_long_homopolymeres_aligned_unsorted_chr1_248745192_248765192.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
#
# find_homopolymers(cram_path, output_path, reference_path)


# cram_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned.bam"
# output_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_long_homopolymeres_aligned_unsorted.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)

# cram_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_chr1_3736056_3936056.bam"
# output_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_long_homopolymeres_aligned_unsorted_chr1_3736056_3936056.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)

# cram_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_chrY_19000000_23000000.bam"
# output_path = "/data/deepvariants/gridss/030945_assembly_ua_realigned_long_homopolymeres_aligned_unsorted_chrY_19000000_23000000.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)






parser = argparse.ArgumentParser(description='Find homopolymeres in reads and realign them using Smith-Waterman algorithm with affine gap penalties and soft clipping')
parser.add_argument('--input', required=True, help='The input CRAM file')
parser.add_argument('--output', required=True, help='The output CRAM file')
parser.add_argument('--reference', required=True, help='The reference genome FASTA file')
parser.add_argument('--homopolymer_length', type=int, default=10, help='The length of homopolymeres to search for')
parser.add_argument("--n_jobs", help="n_jobs of parallel on contigs", type=int, default=-1)
args = parser.parse_args()

MIN_CONTIG_LENGTH = 100000
with pysam.AlignmentFile(args.input, "rc") as cram_file:
    # Get the list of contig names
    contigs = cram_file.references
    # Get the list of contig lengths
    contig_lengths = cram_file.lengths
    contigs = [
        contigs[i] for i in range(len(contigs)) if contig_lengths[i] > MIN_CONTIG_LENGTH
    ]


    Parallel(n_jobs=args.n_jobs, max_nbytes=None)(
        delayed(find_homopolymers)(
            args.input, f"{args.output}{contig}.bam", args.reference, args.homopolymer_length, contig
        )
        for contig in contigs
    )
    # sort and index the contig files
    for contig in contigs:
        pysam.sort("-o", f"{args.output}{contig}_sorted.bam", f"{args.output}{contig}.bam")
        pysam.index(f"{args.output}{contig}_sorted.bam")

    # merge the contig file together
    with pysam.AlignmentFile(args.output, mode='wb', header=cram_file.header) as output:
        for contig in contigs:
            with pysam.AlignmentFile(f"{args.output}{contig}_sorted.bam") as contig_file:
                for read in contig_file:
                    output.write(read)


    # remove the contig files
    for contig in contigs:
        os.remove(f"{args.output}{contig}.bam")
        os.remove(f"{args.output}{contig}_sorted.bam")
        os.remove(f"{args.output}{contig}_sorted.bam.bai")

    pysam.index(args.output)



# find_homopolymers(args.input, args.output, args.reference, args.homopolymer_length)