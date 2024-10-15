import sys
import pysam

def load_record_ids(record_id_file: str) -> set:
    record_ids = set()
    with open(record_id_file) as file:
        for line in file:
            record_ids.add(line.strip())

    return record_ids

def filter_bam_file(bam_file: str, bam_out_file: str, record_ids: set):

    with pysam.AlignmentFile(bam_file) as bam, pysam.AlignmentFile(bam_out_file, "wb", template=bam) as out:
        processed = 0
        included = 0
        for record in bam:
            processed += 1

            if record.query_name in record_ids:
                out.write(record)
                included += 1

        print(f"Included {included} records")

    # Sort and index output bam file
    pysam.sort("-o", bam_out_file, bam_out_file)
    pysam.index(bam_out_file)

if __name__ == "__main__":
    record_ids = load_record_ids(sys.argv[1])
    bam_file = sys.argv[2]
    bam_out_file = sys.argv[3]
    filter_bam_file(bam_file, bam_out_file, record_ids)
