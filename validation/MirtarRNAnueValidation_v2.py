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
    with open(bed_file) as file, open(output_file, "w") as out:

        file.readline()

        for line1, line2 in zip(file, file):
            region1 = get_region_from_bed_line(line1)
            region2 = get_region_from_bed_line(line2)

            if len(list(gff_db.region(region1, completely_within=False, featuretype="miRNA"))) != 0 or len(list(gff_db.region(region2, completely_within=False, featuretype="miRNA"))) != 0:
                features1 = list(gff_db.region(region1, completely_within=False))
                features2 = list(gff_db.region(region2, completely_within=False))

                features1_ids = []
                feature1_miRNA = [feature.attributes["product"] for feature in list(gff_db.region(region1, completely_within=False, featuretype="miRNA"))]

                for feature in features1:
                    if feature.featuretype == "region":
                        continue

                    if feature.featuretype == "gene" and "ID" in feature.attributes:
                        features1_ids.append(feature.attributes["ID"])

                features2_ids = []
                feature2_miRNA = [feature.attributes["product"] for feature in list(gff_db.region(region2, completely_within=False, featuretype="miRNA"))]

                for feature in features2:
                    if feature.featuretype == "region":
                        continue

                    if feature.featuretype == "gene" and "ID" in feature.attributes:
                        features2_ids.append(feature.attributes["ID"])

                out.write(f"Region 1: {region1}, Region 2: {region2}; ")
                out.write(f"miRNA 1: {feature1_miRNA}, miRNA 2: {feature2_miRNA}; ")
                out.write(f"Features 1: {features1_ids}; ")
                out.write(f"Features 2: {features2_ids}\n")




if __name__ == "__main__":
    db_file = sys.argv[4]
    db = load_gff_db(sys.argv[1], db_file)
    bed_file = sys.argv[2]
    output_file = sys.argv[3]
    process_bed_file(bed_file, db, output_file)
