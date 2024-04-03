import pysam
import pyfaidx
import swalign
import io
import re
from Bio import Align

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

def remove_long_homopolymers(sequence, homopolymer_length=10, read_name=""):
    i = 0
    edited_sequence = ""
    start_end_del_tuples = []
    while i < len(sequence):
        count = 1
        while i + count < len(sequence) and sequence[i] == sequence[i + count]:
            count += 1

        if count > homopolymer_length:
            print(
                f"Read {read_name} contains a homopolymer: {sequence[i] * count} at position {i} of the length {count}")
            edited_sequence += sequence[i] * homopolymer_length
            i += count
            start_end_del_tuples.append((i - (count - homopolymer_length) , i))
        elif count > 1:
            edited_sequence += sequence[i] * count
            i += count
        else:
            edited_sequence += sequence[i]
            i += 1
    return edited_sequence, start_end_del_tuples


def find_homopolymers(cram_path, output_path, reference_path, homopolymer_length=10):
    # Open the CRAM file
    reference = pyfaidx.Fasta(reference_path)

    # match = 1
    # mismatch = -4
    # gap_penalty = -6
    # gap_extension_penalty = -1
    # sc_penalty = -5

    #scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    #sw = swalign.LocalAlignment(scoring, gap_penalty=gap_penalty, gap_extension_penalty=gap_extension_penalty)


    # Create a new PairwiseAligner object
    aligner = Align.PairwiseAligner()

    aligner.match_score = 1
    aligner.mismatch_score = -4
    # Set gap scoring
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -1

    # Set end gap scores to 0.0 to not penalize them
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = -5




    with pysam.AlignmentFile(cram_path) as cram:
        with pysam.AlignmentFile(output_path, mode='wb', header=cram.header) as output:
            # Iterate over the reads in the CRAM file

            for read in cram.fetch():
                sequence = read.query_sequence
                if sequence is None:
                    continue  # Skip reads without a query sequence

                # Search for homopolymers in the sequence
                edited_sequence, start_end_del_tuples = remove_long_homopolymers(read.query_sequence, homopolymer_length, read.query_name)

                #Search for homopolymer in the reference
                sc_length = 0
                del_length = 0
                if read.cigartuples[0][0] == 4:
                    # the read is soft clipped
                    sc_length = read.cigartuples[0][1]
                    del_length = start_end_del_tuples[0][1] - start_end_del_tuples[0][0] if len(start_end_del_tuples)>0 else 0
                fa_seq = reference[read.reference_name][read.reference_start - sc_length - del_length:read.reference_start + len(sequence) - sc_length].seq.upper()
                edited_fa_seq, ref_start_end_del_tuples = remove_long_homopolymers(fa_seq, homopolymer_length, "ref of" + read.query_name)


                if len(start_end_del_tuples)>0 or len(ref_start_end_del_tuples)>0:
                    # we found long homopolymer and want to run alignment on that with homopolymere of the length of homopolymer_length

                    print(f"Original sequence:  {sequence}")
                    print(f"Edited sequence:    {edited_sequence}")

                    print(f"Original reference: {fa_seq}")
                    print(f"Edited reference:   {edited_fa_seq}")

                    print("Read cigar: ", read.cigarstring)


                    print("Running alignment on original sequence")
                    score_orig, cigar_orig, start_pos, q_start, r_start = run_alignment_biopyhon(edited_fa_seq, edited_sequence, read.reference_start, sc_length, del_length, aligner)
                    #score_orig, cigar_orig, start_pos = run_alignment(edited_fa_seq, edited_sequence, read.reference_start, sc_length, sw)
                    print(f"score of orig SW:   {score_orig}")

                    # run alignment on reverse complement
                    print("Running alignment on reverse complement sequence")
                    rev_comp_edited_sequence = reverse_complement(edited_sequence)
                    score_rev, cigar_rev, start_pos_rev, q_start_rev, r_start_rev = run_alignment_biopyhon(edited_fa_seq, edited_sequence, read.reference_start, sc_length, del_length, aligner)
                    #score_rev, cigar_rev, start_pos_rev = run_alignment(edited_fa_seq, rev_comp_edited_sequence,read.reference_start, 0, sw)#TODO

                    print(f"score of rev SW:    {score_rev}")

                    if score_orig >= score_rev:
                        print("Original sequence is better")
                        score = score_orig
                        cigar = cigar_orig
                        updated_seq = edited_sequence
                    else:
                        print("Reverse complement sequence is better")
                        score = score_rev
                        cigar = cigar_rev
                        updated_seq = rev_comp_edited_sequence
                        start_end_del_tuples = [(len(sequence) - end, len(sequence) - start) for start, end in start_end_del_tuples]
                        ref_start_end_del_tuples = [(len(fa_seq) - end, len(fa_seq) - start) for start, end in ref_start_end_del_tuples]

                    print(f"Final score:                      {score}")

                    print("original start pos: ", read.reference_start)
                    print("updated start pos: ", start_pos)
                    print("original seq length: ", len(sequence))
                    print("Read cigar (length): ", calculate_sequence_length_by_cigar(read.cigarstring))
                    print("updated seq length: ", len(updated_seq))

                    # Insert deletions into the CIGAR string
                    updated_cigar = cigar
                    print(f"Final cigar before change:        {cigar}")
                    print("Final cigar before change (length): ", calculate_sequence_length_by_cigar(cigar))
                    print(f"seq start_end_del_tuples: {start_end_del_tuples}")
                    print(f"                             {updated_cigar}")
                    for start_del, end_del in start_end_del_tuples:
                        updated_cigar = insert_operation_into_cigar(updated_cigar, start_del, end_del - start_del, 'I')
                        print(f"start_del: {start_del}, end_del: {end_del}, cigar_length: {calculate_sequence_length_by_cigar(updated_cigar)} ", updated_cigar)

                    updated_cigar_D = updated_cigar

                    print(f"ref start_end_del_tuples: {ref_start_end_del_tuples}")
                    print(f"                             {updated_cigar}")
                    for start_del, end_del in ref_start_end_del_tuples:
                        updated_cigar = insert_operation_into_cigar(updated_cigar, start_del - r_start + q_start, end_del - start_del, 'D')
                        print(f"start_del: {start_del}, end_del: {end_del}, cigar_length: {calculate_sequence_length_by_cigar(updated_cigar)} ", updated_cigar)





                    print(f"updated_cigar_D                   {updated_cigar_D}")
                    print("Final cigar before change_D (length): ", calculate_sequence_length_by_cigar(updated_cigar_D))
                    print(f"Final cigar after change:         {updated_cigar}")
                    print("Final cigar after change (length): ", calculate_sequence_length_by_cigar(updated_cigar))
                    print("cigartuples: ", cigar_string_to_cigartuples(updated_cigar))

                    read.cigar = cigar_string_to_cigartuples(updated_cigar)
                    read.reference_start = start_pos

                    assert(len(sequence) == calculate_sequence_length_by_cigar(updated_cigar))


                    output.write(read)

                    print("#########")
                    print("\n")
                else:
                    print(f"Original sequence:  {sequence}")
                    print("No homopolymer found")
                    print("#########")
                    print("\n")
                    output.write(read)


def pad_alignment_with_sc(cigar, q_pos, q_end, r_pos, r_end, expected_length, start_pos, sc_length):
    """Pad the CIGAR string with soft-clipping operations to match the expected length."""
    # Regular expression to find numbers followed by CIGAR operation characters
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')

    # Find all matches of the pattern in the CIGAR string
    matches = pattern.findall(cigar)

    # Sum up the lengths of operations that consume the query sequence
    consumed_length = 0
    for count, op in matches:
        if op in "MIS=X":
            consumed_length += int(count)

    # Calculate the length of soft-clipping to add
    if q_pos < sc_length:
        # Means we should start ref a bit earlier
        start_pos = start_pos - sc_length + r_pos
        #start_pos = start_pos - (sc_length + (0 if r_pos == 0 else (r_pos -1)))
        #start_pos = start_pos - (sc_length + (0 if r_pos == 0 else (r_pos - 1)))



    padded_cigar = cigar
    if q_pos > 0:
        # Add the soft-clipping to the beginning of the CIGAR string
        padded_cigar = f"{q_pos}S{padded_cigar}"

    if q_end < expected_length: # TODO
        # Add the soft-clipping to the end of the CIGAR string
        padded_cigar = f"{padded_cigar}{expected_length-q_end}S"

    return padded_cigar, start_pos


def convert_aligned_to_cigar(aligned, seq1_len, seq2_len):
    target_aligned, query_aligned = aligned
    cigar = []
    adjusted_start_pos = 0  # This will store the adjusted start position based on initial target insertions

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

        # Adjust for gaps at the beginning of alignment segments
        elif t_gap > 0 or q_gap > 0:
            # If there's a gap in target (deletion relative to query)
            if t_gap > q_gap:
                add_op('D', t_gap - q_gap)
            # If there's a gap in query (insertion relative to target)
            elif q_gap > t_gap:
                add_op('I', q_gap - t_gap)

        # Add matched segment
        segment_length = min(t_end - t_start, q_end - q_start)
        add_op('M', segment_length)

        # Update last processed positions
        last_target_end, last_query_end = t_end, q_end

    # Handle gaps after the last aligned segment
    if last_target_end < seq1_len:
        add_op('D', seq1_len - last_target_end)
    if last_query_end < seq2_len:
        add_op('S', seq2_len - last_query_end)  # Treat remaining query as soft clipped

    print(f"adjusted_start_pos: {adjusted_start_pos}")
    return adjusted_start_pos, ''.join(cigar), query_aligned[0][0], target_aligned[0][0]

def run_alignment_biopyhon(fa_seq, sequence, start_pos, sc_length, del_length, aligner):
    # Perform the alignment between two sequences
    for alignment in aligner.align(fa_seq, sequence ):
        # Print each alignment's score and the alignment itself
        print("Score = %.1f:" % alignment.score)
        start_pos_adjust, cigar, q_start, r_start = convert_aligned_to_cigar(alignment.aligned, len(fa_seq), len(sequence))
        print(F"aliged = {alignment.aligned}")
        print(f"Cigar = {cigar}")
        print(alignment)
        start_pos = start_pos - sc_length - del_length + start_pos_adjust
        return alignment.score, cigar, start_pos, q_start, r_start


def run_alignment(fa_seq, sequence, start_pos, sc_length, sw):

                alignment = sw.align(fa_seq, sequence)
                with io.StringIO() as file:
                    alignment.dump(out=file)
                    out = file.getvalue()
                    match = re.search(r'CIGAR: ([\dMIDNSHP=X]+)', out)
                    cigar = match.group(1)
                    match = re.search(r'Score: (\d+)', out)
                    score = int(match.group(1))
                    print(out)
                    print(f"alignment.q_pos: {alignment.q_pos} alignment.q_end: {alignment.q_end}")
                    print(f"alignment.r_pos: {alignment.r_pos} alignment.r_end: {alignment.r_end}")
                    padded_cigar, start_pos = pad_alignment_with_sc(cigar, alignment.q_pos, alignment.q_end, alignment.r_pos, alignment.r_end, len(sequence), start_pos, sc_length)
                    return score, padded_cigar, start_pos
def reverse_complement(seq):
    """Return the reverse complement of the given DNA sequence."""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    # Generate the complement sequence
    complement_seq = ''.join(complement_dict[nuc] for nuc in seq)
    # Reverse the complement sequence
    reverse_comp_seq = complement_seq[::-1]
    return reverse_comp_seq



def insert_operation_into_cigar(cigar, position, op_size, op_type):
    assert op_type in ['D', 'I'], "op_type must be 'D' for deletion or 'I' for insertion"

    # Split the CIGAR string into its components
    cigar_components = re.findall(r'\d+[MIDNSHP=X]', cigar)

    # Convert components into a list of tuples (count, operation)
    parsed_cigar = [(int(comp[:-1]), comp[-1]) for comp in cigar_components]

    new_cigar_ops = []
    accumulated_length = 0
    operation_inserted = False
    i=0
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
            if (not skip_operation) and (accumulated_length < position <= accumulated_length + count):
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

    return updated_cigar

# # usage
# cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9_69449415_69450719.bam"
# output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9_69449415_69450719.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)
#

cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9.bam"
output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9.bam"
reference_path = "/data/Homo_sapiens_assembly38.fasta"
find_homopolymers(cram_path, output_path, reference_path)


# cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9_275409_339530.bam"
# output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9_275409_339530.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)

# cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9_1_275409.bam"
# output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9_1_275409.bam"
# reference_path = "/data/Homo_sapiens_assembly38.fasta"
# find_homopolymers(cram_path, output_path, reference_path)




# with pysam.AlignmentFile("/data/deepvariants/gridss/030945_merged_assembly_chr9_69449415_69450719.bam") as infile:
#     with pysam.AlignmentFile('/data/deepvariants/gridss/test.bam', mode='w', header=infile.header) as outfile:
#          for rec in infile:
#             #rec.cigarstring= '115M1I1D'
#             rec.query_sequence = 'AGAAAGACACTTCCCTGGGGACACCAACCCAGATGAGTTCCTGTCTTCTCAGCATTCCGCATATTTGGAGTTTTTAAGAAATGAATTCACACAGGTCTACACTCTTTTGTAATTCTCTCGTTTCACATAAGCAAACTTGCCTCAGCACACAACCATGAGGACCACCAGTTTTTTTTTTATTCATTTCGTCC'
#             rec.cigar = [(0, 115), (1, 1), (2, 1), (0, 62), (2, 10), (1, 1), (0, 2), (1, 1),(2, 1),(1, 1),(0, 3),(4, 5)]
#             #rec.cigar = [(0, 115), (1, 1), (2, 1), (0, 62), (2, 10), (1, 1), (0, 2), (1, 1), (2, 1), (1, 1),(0, 3),(0, 5)]
#             #rec.query_sequence = 'A'*100
#             outfile.write(rec)
