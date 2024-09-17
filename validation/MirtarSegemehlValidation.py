import sys
import subprocess
from pathlib import Path


def filter_sam_for_split_reads(sam_file: Path, out_bam_file: Path):
    command: list[str] = [
        "samtools",
        "view",
        "-bh",
        "-o",
        out_bam_file.as_posix(),
        "-@",
        "8",
        sam_file.as_posix(),
    ]

    return_value = subprocess.run(command, check=True)

    if return_value.returncode != 0:
        print(f"Error: {return_value.stderr}")
        sys.exit(1)


def intersect_bam_file_with_mirtar_annotation(
    bam_file: Path,
    mirtar_annotation_file: Path,
    out_bam_overlaps_file: Path,
    out_overlaps_count_file: Path,
):
    command_intersect: list[str] = [
        "bedtools",
        "intersect",
        "-split",
        "-sortout",
        "-u",
        "-a",
        bam_file.as_posix(),
        "-b",
        mirtar_annotation_file.as_posix(),
    ]

    out_bam_overlaps = open(out_bam_overlaps_file, "w")
    return_value = subprocess.run(
        command_intersect, check=True, stdout=out_bam_overlaps
    )
    out_bam_overlaps.close()

    if return_value.returncode != 0:
        print(f"Error: {return_value.stderr}")
        sys.exit(1)

    command_count: list[str] = [
        "bedtools",
        "intersect",
        "-split",
        "-sortout",
        "-c",
        "-a",
        mirtar_annotation_file.as_posix(),
        "-b",
        bam_file.as_posix(),
    ]

    out_overlaps_count = open(out_overlaps_count_file, "w")
    return_value = subprocess.run(command_count, check=True, stdout=out_overlaps_count)
    out_overlaps_count.close()

    if return_value.returncode != 0:
        print(f"Error: {return_value.stderr}")
        sys.exit(1)


def main():
    sam_file: Path = Path(sys.argv[1])
    mirtar_annotation_file: Path = Path(sys.argv[2])
    out_dir: Path = Path(sys.argv[3])

    out_dir.mkdir(parents=True, exist_ok=True)

    sample_name: str = sam_file.stem

    out_bam_file: Path = out_dir / f"{sample_name}_splits_filtered.bam"

    filter_sam_for_split_reads(sam_file, out_bam_file)

    out_bam_overlaps_file: Path = (
        out_dir / f"{sample_name}_splits_filtered_overlaps.bam"
    )
    out_overlaps_count_file: Path = (
        out_dir / f"{sample_name}_splits_filtered_overlaps_count.tsv"
    )

    intersect_bam_file_with_mirtar_annotation(
        out_bam_file,
        mirtar_annotation_file,
        out_bam_overlaps_file,
        out_overlaps_count_file,
    )


if __name__ == "__main__":
    main()
