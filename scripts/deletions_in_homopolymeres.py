import pysam
import pyfaidx
import swalign
import io
import re

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

def find_homopolymers(cram_path, output_path, reference_path, homopolymer_length=10):
    # Open the CRAM file
    reference = pyfaidx.Fasta(reference_path)

    match = 1
    mismatch = -4
    gap_penalty = -1
    gap_extension_penalty = -6


    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=gap_penalty, gap_extension_penalty=gap_extension_penalty)


    with pysam.AlignmentFile(cram_path, "rc", reference_filename=reference_path) as cram:
        with pysam.AlignmentFile(output_path, "wc", template=cram) as output:
            # Iterate over the reads in the CRAM file

            for read in cram.fetch():
                sequence = read.query_sequence
                if sequence is None:
                    continue  # Skip reads without a query sequence

                # Search for homopolymers in the sequence
                i = 0
                start_del=-1
                end_del=-1
                found = False
                edited_sequence = ""
                while i < len(sequence):
                    count = 1
                    while i + count < len(sequence) and sequence[i] == sequence[i + count] and found == False:
                        count += 1

                    if count >= homopolymer_length:
                        print(f"Read {read.query_name} contains a homopolymer: {sequence[i] * count} at position {i} of the length {count}")
                        edited_sequence += sequence[i] * homopolymer_length
                        i += count
                        found = True
                        start_del = i-count+homopolymer_length
                        end_del = i
                    elif count > 1:
                        edited_sequence += sequence[i] * count
                        i+=count
                    else:
                        edited_sequence += sequence[i]
                        i+=1

                if found: # we found long homopolymer and want to run alignment on that with homopolymere of the length of homopolymer_length
                    fa_seq = reference[read.reference_name][read.reference_start-1:read.reference_start-1+len(sequence)].seq.upper()
                    edited_fa_seq = fa_seq[:start_del] + fa_seq[end_del:]
                    print(f"Original sequence:  {sequence}")
                    print(f"Edited sequence:    {edited_sequence}")

                    print(f"Original reference: {fa_seq}")
                    print(f"Edited reference:   {edited_fa_seq}")

                    print("Read cigar: ", read.cigarstring)



                    score_orig, cigar_orig = run_alignment(edited_fa_seq, edited_sequence, sw)
                    print(f"score of orig SW:   {score_orig}")

                    # run alignment on reverse complement
                    rev_comp_edited_fa_seq = reverse_complement(edited_fa_seq)
                    rev_comp_edited_sequence = reverse_complement(edited_sequence)
                    score_rev, cigar_rev = run_alignment(rev_comp_edited_fa_seq, rev_comp_edited_sequence, sw)

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
                        start_del = len(sequence) - end_del
                        end_del = len(sequence) - start_del

                    print(f"Final score:                      {score}")


                    updated_cigar = insert_deletion_into_cigar(cigar, start_del, end_del - start_del)

                    print("original seq length: ", len(sequence))
                    print("Read cigar (length): ", calculate_sequence_length_by_cigar(read.cigarstring))
                    print("updated seq length: ", len(updated_seq))
                    print(f"Final cigar before change:        {cigar}")
                    print("Final cigar before change (length): ", calculate_sequence_length_by_cigar(cigar))
                    print(f"Final cigar after change:         {updated_cigar}")
                    print("Final cigar after change (length): ", calculate_sequence_length_by_cigar(updated_cigar))
                    print("cigartuples: ", cigar_string_to_cigartuples(updated_cigar))

                    read.cigarstring = updated_cigar
                    read.query_sequence = sequence


                    output.write(read)

                    print("#########")
                    print("\n")
                else:
                    print(f"Original sequence:  {sequence}")
                    print("No homopolymer found")
                    print("#########")
                    print("\n")
                    output.write(read)



def run_alignment(fa_seq, sequence, sw):
                alignment = sw.align(fa_seq, sequence)
                with io.StringIO() as file:
                    alignment.dump(out=file)
                    out = file.getvalue()
                    match = re.search(r'CIGAR: ([\dMIDNSHP=X]+)', out)
                    cigar = match.group(1)
                    match = re.search(r'Score: (\d+)', out)
                    score = int(match.group(1))
                    print(out)
                    return score, cigar
def reverse_complement(seq):
    """Return the reverse complement of the given DNA sequence."""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    # Generate the complement sequence
    complement_seq = ''.join(complement_dict[nuc] for nuc in seq)
    # Reverse the complement sequence
    reverse_comp_seq = complement_seq[::-1]
    return reverse_comp_seq


def merge_operations(ops):
    """Merge consecutive deletions and adjust for insertion followed by deletion."""
    merged_ops = []
    skip_next = False

    for i, (count, op) in enumerate(ops[:-1]):
        if skip_next:
            skip_next = False
            continue

        # Look ahead to next operation
        next_count, next_op = ops[i + 1]

        # Merge consecutive deletions
        if op == 'D' and next_op == 'D':
            merged_ops.append((count + next_count, 'D'))
            skip_next = True
        # Adjust for insertion followed by deletion
        elif op == 'I' and next_op == 'D':
            if count > next_count:
                merged_ops.append((count - next_count, 'I'))
            elif count < next_count:
                merged_ops.append((next_count - count, 'D'))
            # If equal, they cancel each other, no operation is appended
            skip_next = True
        else:
            merged_ops.append((count, op))

    if not skip_next:  # Ensure the last operation is added if not already handled
        merged_ops.append(ops[-1])

    return merged_ops


def insert_deletion_into_cigar(cigar, position, deletion_size):
    import re

    # Split the CIGAR string into its components
    cigar_components = re.findall(r'\d+[MIDNSHP=X]', cigar)

    # Convert components into a list of tuples (count, operation)
    parsed_cigar = [(int(comp[:-1]), comp[-1]) for comp in cigar_components]

    new_cigar_ops = []
    accumulated_length = 0
    deletion_inserted = False

    for count, op in parsed_cigar:
        # Only proceed if we haven't inserted the deletion yet
        if not deletion_inserted:
            # Operations that consume the reference (for position tracking)
            if op in 'MDN=X':
                # If the position falls within this operation's span
                if accumulated_length < position <= accumulated_length + count:
                    # Calculate the split before and after the deletion
                    before_length = position - accumulated_length
                    after_length = count - before_length

                    # Add the operation before the deletion, if any
                    if before_length > 0:
                        new_cigar_ops.append((before_length, op))

                    # Insert the deletion
                    new_cigar_ops.append((deletion_size, 'D'))
                    deletion_inserted = True

                    # Add what's left of the operation after the deletion, if any
                    if after_length > 0:
                        new_cigar_ops.append((after_length, op))
                else:
                    new_cigar_ops.append((count, op))
                accumulated_length += count
            else:  # For operations that do not consume the reference
                new_cigar_ops.append((count, op))
        else:
            # Add all remaining operations after the deletion has been inserted
            new_cigar_ops.append((count, op))

    # Convert back to CIGAR string format
    updated_cigar = ''.join(f"{count}{op}" for count, op in new_cigar_ops)

    return updated_cigar

# usage
cram_path = "/data/deepvariants/gridss/030945_merged_assembly_chr9_69449415_69450719.bam"
output_path = "/data/deepvariants/gridss/homopolymeres_030945_merged_assembly_chr9_69449415_69450719.bam"
reference_path = "/data/Homo_sapiens_assembly38.fasta"
find_homopolymers(cram_path, output_path, reference_path)