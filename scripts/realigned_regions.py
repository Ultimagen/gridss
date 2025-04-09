import pysam
from collections import defaultdict
import sys
import argparse

logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__ if __name__ != "__main__" else "realigned_regions")

def process_bam_files(before_bam, after_bam, output_bed):
    # Open BAM files
    before = pysam.AlignmentFile(before_bam, "rb")
    after = pysam.AlignmentFile(after_bam, "rb")

    # Dictionary to store changed regions
    changed_regions = defaultdict(list)

    def add_region(chrom, start, end):
        changed_regions[chrom].append((start, end))

    # Read all entries from "before" BAM into dictionary
    before_dict = defaultdict(list)
    for read in before:
        key = (read.query_name, read.flag)
        before_dict[key].append((read.reference_name, read.reference_start, read.reference_end))

    # Process "after" BAM file
    for read in after:
        key = (read.query_name, read.flag)
        current_position = (read.reference_name, read.reference_start, read.reference_end, read.cigartuples)

        if key in before_dict:
            previous_positions = before_dict[key][0]

            if (current_position[0] != previous_positions[0] or abs(current_position[1] - previous_positions[1]) > 500 or \
                    abs(current_position[2] - previous_positions[2]) > 500) and \
                    (any(op in {1, 2, 4} and length > 20 for op, length in (read.cigartuples or [])) or
                    any(op in {1, 2, 4} and length > 20 for op, length in (previous_positions[3] or []))):
                # Read changed location
                add_region(previous_positions[0], previous_positions[1], previous_positions[2])
                add_region(current_position[0], current_position[1], current_position[2])

            del before_dict[key]  # Remove matched read
        else:
            # New read
            add_region(read.reference_name, read.reference_start, read.reference_end)

    # Remaining reads in "before_dict" are missing in "after.bam"
    for (query_name, flag), positions in before_dict.items():
        for ref_name, ref_start, ref_end, _ in positions:
            add_region(ref_name, ref_start, ref_end)

    # Merge overlapping regions and sort
    merged_regions = []
    for chrom in changed_regions:
        intervals = sorted(changed_regions[chrom])
        merged = []
        start, end = intervals[0]

        for i in range(1, len(intervals)):
            s, e = intervals[i]
            if s <= end:  # Overlapping regions
                end = max(end, e)
            else:
                merged.append((start, end))
                start, end = s, e
        merged.append((start, end))
        merged_regions.extend([(chrom, s, e) for s, e in merged])

    # Calculate total length of changed regions
    total_length = sum(e - s for _, s, e in merged_regions)

    # Write sorted merged regions to BED file
    with open(output_bed, "w") as bed_file:
        for chrom, start, end in sorted(merged_regions):
            bed_file.write(f"{chrom}\t{start}\t{end}\n")

    # Close BAM files
    before.close()
    after.close()

    logger.info(f"Merged and ordered regions saved in {output_bed}.")
    logger.info(f"Total changed region length: {total_length} bases.")


parser = argparse.ArgumentParser(description="Extract realigned regions from BAM files")
parser.add_argument("--before_bam", help="Original BAM file")
parser.add_argument("--after_bam", help="Realigned BAM file")
parser.add_argument("--output_bed", help="Output BED file")
args = parser.parse_args()

process_bam_files(args.before_bam, args.after_bam, args.output_bed)