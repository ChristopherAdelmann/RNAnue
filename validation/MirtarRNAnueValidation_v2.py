# 1. Load gff feature file which is filtered by mrtardb ids
# 2. Iterate over bed file in pairs of 2 lines and check if at least one of the regions overlaps with the gff feature
# 3. If overlap is found extract both features and write them to the output file

import sys
import gffutils
import os.path

def load_gff_db(gff_file: str, db_filename: str) -> gffutils.FeatureDB:
    print(f"Creating database for {gff_file}")

    if os.path.isfile(db_filename):
        print("Database already exists")
        return gffutils.FeatureDB(db_filename)

    db = gffutils.create_db(gff_file, dbfn=db_filename, merge_strategy="merge")
    print("Database created")
    return db

def get_region_from_bed_line(line: str):
    bed_fields = line.strip().split("\t")
    return (bed_fields[0], int(bed_fields[1]), int(bed_fields[2]))

def process_bed_file(bed_file: str, gff_db: gffutils.FeatureDB, output_file: str):
    print(f"Processing bed file {bed_file}")
    with open(bed_file) as file:
        for line1, line2 in zip(file, file):
            region1 = get_region_from_bed_line(line1)
            region2 = get_region_from_bed_line(line2)

            for feature in gff_db.region(region1, completely_within=False, featuretype="miRNA"):
                print(f"Found feature {feature.id} overlapping with {region1}")


            for feature in gff_db.region(region2, completely_within=False, featuretype="miRNA"):
                print(f"Found feature {feature.id} overlapping with {region2}")

if __name__ == "__main__":
    db_file = sys.argv[4]
    db = load_gff_db(sys.argv[1], db_file)
    bed_file = sys.argv[2]
    output_file = sys.argv[3]
    process_bed_file(bed_file, db, output_file)
