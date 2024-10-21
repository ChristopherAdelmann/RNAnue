import sys
import pysam
from pathlib import Path
import subprocess

bam_file = Path(sys.argv[1])
mirtar_gff_file = Path(sys.argv[2])
db_file = Path(sys.argv[3])
out_dir = Path(sys.argv[4])
mirtar_file = Path(sys.argv[5])

# Extract split alignments from raw alignment file (filter for XJ tag)

print("Extracting split alignments from raw alignment file")

split_bam_file = out_dir / "split.bam"
command = [
    "samtools",
    "view",
    "-hb",
    "-d", "XJ:2",
    "-@", "5",
    "-o", split_bam_file.as_posix(),
    bam_file.as_posix()]

return_value = subprocess.run(command, check=True)

if return_value.returncode != 0:
    print(f"Error: {return_value.stderr}")
    sys.exit(1)

# Intersect miRNA gff with alignment file (bedtools intersect)

print("Intersecting miRNA gff with alignment file")

intersected_file = out_dir / "intersected.bam"
command_intersect = [
    "bedtools",
    "intersect",
    "-split",
    "-abam", split_bam_file.as_posix(),
    "-b", mirtar_gff_file.as_posix(),
    "-u"]

out_bam_overlaps = open(intersected_file, "w")
return_value = subprocess.run(
    command_intersect, check=True, stdout=out_bam_overlaps, stderr=subprocess.PIPE)
out_bam_overlaps.close()

if return_value.returncode != 0:
    print(f"Error: {return_value.stderr}")
    sys.exit(1)

# Extract read IDs from intersected alignments

print("Extracting read IDs from intersected alignments")

read_ids = set()
with pysam.AlignmentFile(intersected_file.as_posix()) as bam:
    for record in bam:
        read_ids.add(record.query_name)

# Filter split alignment file by read IDs

print("Filtering split alignment file by read IDs")

filtered_bam_file = out_dir / "filtered.bam"
with pysam.AlignmentFile(split_bam_file.as_posix()) as bam, pysam.AlignmentFile(filtered_bam_file.as_posix(), "wb", template=bam) as out:
    for record in bam:
        if record.query_name in read_ids:
            out.write(record)

# Sort filtered file by read ID

print("Sorting filtered file by read ID")

sorted_bam_file = out_dir / "sorted_filterd.bam"

command = [
    "samtools",
    "sort",
    "-n",
    "-@",
    "5",
    "-o",
    sorted_bam_file.as_posix(),
    filtered_bam_file.as_posix()]

return_value = subprocess.run(command, check=True)

if return_value.returncode != 0:
    print(f"Error: {return_value.stderr}")
    sys.exit(1)

# Convert sorted file to bed format

print("Converting sorted file to bed format")

bed_file = out_dir / "sorted_filtered.bed"

command = [
    "bedtools",
    "bamtobed",
    "-split",
    "-i", sorted_bam_file.as_posix()
]

out_bed = open(bed_file, "w")
return_value = subprocess.run(
    command, check=True, stdout=out_bed)
out_bed.close()

if return_value.returncode != 0:
    print(f"Error: {return_value.stderr}")
    sys.exit(1)

# Run MirtarRNAnueValidation_v2.py with the bed file (relative to the script in the same directory)

print("Running MirtarRNAValidation_v2.py with bed file")

log_file = out_dir / "log_mirtar_validation.txt"

current_dir = Path(__file__).parent
mirtar_validation_script = current_dir / "MirtarRNAnueValidation_v2.py"

command = [
    "python3",
    mirtar_validation_script.as_posix(),
    mirtar_gff_file.as_posix(),
    bed_file.as_posix(),
    log_file.as_posix(),
    db_file.as_posix(),
    mirtar_file.as_posix()
]

return_value = subprocess.run(command, check=True)

if return_value.returncode != 0:
    print(f"Error: {return_value.stderr}")
    sys.exit(1)
