#!/usr/bin/env python3
# Description: This script processes two BAM files to identify realigned regions and outputs them in BED format.

import pysam
from collections import defaultdict
import argparse
import logging



POSITION_DELTA_THRESHOLD = 500

def parse_args():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Extract realigned regions from BAM files")
    parser.add_argument("--before_bam", required=True, type=str, help="Original BAM file before realignment")
    parser.add_argument("--after_bam", required=True, type=str, help="BAM file after realignment")
    parser.add_argument("--output_bed", required=True, type=str, help="Output BED file for realigned regions")
    args = parser.parse_args()
    return args


def run():
    args = parse_args()
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    logger = logging.getLogger(__name__ if __name__ != "__main__" else "realigned_regions")
    before_bam = args.before_bam
    after_bam = args.after_bam
    output_bed = args.output_bed
    # Open BAM files
    with pysam.AlignmentFile(before_bam, "rb") as before:
        with pysam.AlignmentFile(after_bam, "rb") as after:
            logger.info(f"Finding realigned regions between {before_bam} and {after_bam}")
            # Dictionary to store changed regions
            changed_regions = defaultdict(list)

            def add_region(chrom, start, end):
                changed_regions[chrom].append((start, end))

            # Read all entries from "before" BAM into dictionary
            before_dict = defaultdict(list)
            for read in before:
                key = (read.query_name, read.flag)
                before_dict[key].append((read.reference_name, read.reference_start, read.reference_end, read.cigartuples))

            # Process "after" BAM file
            for read in after:
                key = (read.query_name, read.flag)
                current_position = (read.reference_name, read.reference_start, read.reference_end, read.cigartuples)

                if key in before_dict:
                    previous_positions = before_dict[key][0]

                    if (current_position[0] != previous_positions[0] or
                            abs(current_position[1] - previous_positions[1]) > POSITION_DELTA_THRESHOLD or
                            abs(current_position[2] - previous_positions[2]) > POSITION_DELTA_THRESHOLD):
                        # Read changed location
                        add_region(previous_positions[0], previous_positions[1], previous_positions[2])
                        logger.debug(f"Read {read.query_name} changed location from {previous_positions} to {current_position}")
                        add_region(current_position[0], current_position[1], current_position[2])

                    del before_dict[key]  # Remove matched read
                else:
                    # New read
                    add_region(read.reference_name, read.reference_start, read.reference_end)
                    logger.debug(f"New read {read.query_name} at {current_position}")

            # Remaining reads in "before_dict" are missing in "after.bam"
            for (query_name, flag), positions in before_dict.items():
                for ref_name, ref_start, ref_end, _ in positions:
                    add_region(ref_name, ref_start, ref_end)
                    logger.debug(f"Read {query_name} missing in after.bam at {ref_name}:{ref_start}-{ref_end}")

    logger.info(f"Found {len(changed_regions)} changed regions across all chromosomes.")
    logger.info("Merging overlapping regions and sorting them")
    # Merge overlapping regions and sort
    merged_regions = []
    for chrom in changed_regions:
        intervals = sorted(changed_regions[chrom])
        if not intervals:  # Skip if no intervals for this chromosome
            continue
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

    logger.info(f"Merged and ordered regions saved in {output_bed}.")
    logger.info(f"Total changed region length: {total_length} bases.")




if __name__ == "__main__":
    run()