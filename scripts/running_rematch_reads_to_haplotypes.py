import subprocess
import multiprocessing
import os

# Load BED file regions
bed_file = "/data/Runs/bioin_1984_use_both_sides/fp_tp_regions_sorted.bed"
with open(bed_file, "r") as f:
    regions = [line.strip().split() for line in f]


# Define function to run each script instance
def run_script(region):
    chrom, start, end = region
    region_str = f"{chrom}:{start}-{end}"
    output_bam = f"/data/Runs/bioin_1984_use_both_sides/output_region_{region_str}.bam"
    # in case output bam exists - skip
    if os.path.exists(f"{output_bam}.bai"):
        return


    cmd = [
        "/home/ubuntu/miniconda3/envs/genomics.py3/bin/python", "/home/ubuntu/proj/VariantCalling/ugvc/rematch_reads_to_haplotypes.py",
        "--assembly",
        "/data/Runs/bioin_1984_use_both_sides/sv_somatic_colo_assembly_ua_realigned_long_homopolymers_aligned.bam",
        "--reference", "/data/Runs/Homo_sapiens_assembly38.fasta",
        "--output", output_bam,
        "--tumor_crams", "/data/Runs/bioin_1984_use_both_sides/tumor.031865-Lb_2211-Z0048-CTGCCAGACTGTGAT.cram",
        "--germline_crams", "/data/Runs/bioin_1984_use_both_sides/normal.031865-Lb_2212-Z0134-CTTCAGCATACAGAT.cram",
        "--region", region_str
    ]
    print(region_str)
    subprocess.run(cmd)


# Run in parallel with 48 processes
if __name__ == "__main__":
    with multiprocessing.Pool(2) as pool:
        pool.map(run_script, regions)
        pool.close()
        pool.join()


# samtools merge -u - input1.bam input2.bam input3.bam \
#   | samtools view -h - \
#   | sed 's/\tRG:Z:[^\t]*//g' \
#   | awk '!seen[$0]++' \
#   | samtools view -b -o merged_no_RG_unique.bam

# Merge, strip RG tags, deduplicate, and output cleaned BAM
output_bams = [f"/data/Runs/bioin_1984_use_both_sides/output_region_{region[0]}:{region[1]}-{region[2]}.bam" for region in regions]
merged_clean_bam = "/data/Runs/bioin_1984_use_both_sides/merged_output_no_RG_unique.bam"
sorted_bam = "/data/Runs/bioin_1984_use_both_sides/merged_output.bam"

# Merge the BAM files
merge_command = ["samtools", "merge", "-u", "-f", "-", *output_bams]
process = subprocess.Popen(merge_command, stdout=subprocess.PIPE)

# Remove RG tags, deduplicate, and write to final BAM
view_command = ["samtools", "view", "-h", "-"]
sed_command = ["sed", r's/\tRG:Z:[^\t]*//g']
awk_command = ["awk", '!seen[$0]++']
convert_command = ["samtools", "view", "-b", "-o", merged_clean_bam]

# Chain the commands together
p1 = subprocess.Popen(view_command, stdin=process.stdout, stdout=subprocess.PIPE)
p2 = subprocess.Popen(sed_command, stdin=p1.stdout, stdout=subprocess.PIPE, shell=False)
p3 = subprocess.Popen(awk_command, stdin=p2.stdout, stdout=subprocess.PIPE, shell=False)
p4 = subprocess.Popen(convert_command, stdin=p3.stdout)

# Wait for the process to finish
p4.communicate()

# Sort and index the final BAM
subprocess.run(["samtools", "sort", "-o", sorted_bam, merged_clean_bam])
subprocess.run(["samtools", "index", sorted_bam])

# Clean up intermediate files only after confirming they're created
if os.path.exists(merged_clean_bam):
    os.remove(merged_clean_bam)
if os.path.exists(f"{merged_clean_bam}.bai"):
    os.remove(f"{merged_clean_bam}.bai")
