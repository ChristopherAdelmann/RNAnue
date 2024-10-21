# 1. Load gff feature file which is filtered by mrtardb ids
# 2. Iterate over bed file in pairs of 2 lines and check if at least one of the regions overlaps with the gff feature
# 3. If overlap is found extract both features and write them to the output file

from collections import defaultdict
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

def get_mirna_mrna_pair(gff_db: gffutils.FeatureDB, miRNA_region: tuple[str, int, int], mRNA_region: tuple[str, int, int]) -> tuple[list[gffutils.Feature], list[gffutils.Feature]] | None:
    miRNA_features = list(gff_db.region(miRNA_region, completely_within=False, featuretype="miRNA"))
    if len(miRNA_features) == 0:
        return None

    mRNA_features = list(gff_db.region(mRNA_region, completely_within=False, featuretype="gene"))

    return miRNA_features, mRNA_features

def load_mirtar_interactions(mirtar_file: str) -> dict[str, set[str]]:
    interactions: dict[str, set[str]] = defaultdict(set)

    with open(mirtar_file) as file:
        file.readline()

        for line in file:
            fields = line.strip().split("\t")
            interactions[fields[1]].add(fields[3])

    return interactions

def process_bed_file(bed_file: str, gff_db: gffutils.FeatureDB, output_file: str):
    print(f"Processing bed file {bed_file}")
    with open(bed_file) as file, open(output_file, "w") as out:

        interaction_counts: dict[tuple[str, str], int] = defaultdict(int)

        for line1, line2 in zip(file, file):
            read_id1 = line1.strip().split("\t")[3]
            read_id2 = line2.strip().split("\t")[3]

            if read_id1 != read_id2:
                print(f"Read ids do not match: {read_id1} != {read_id2}")

                exit(0)

                continue

            region1 = get_region_from_bed_line(line1)
            region2 = get_region_from_bed_line(line2)

            result1 = get_mirna_mrna_pair(gff_db, region1, region2)

            if result1:
                miRNA, mRNA = result1
                features1_ids = ", ".join([feature.attributes["product"][0] for feature in miRNA])
                features2_ids = ", ".join([feature.attributes["gene"][0] for feature in mRNA])

                out.write(f"Region 1: {region1}, Region 2: {region2}; ")
                out.write(f"miRNA: {features1_ids}, mRNA: {features2_ids}\n")

                for miRNA_feature in miRNA:
                    for mRNA_feature in mRNA:
                        interaction_counts[(miRNA_feature.attributes["product"][0], mRNA_feature.attributes["gene"][0])] += 1

            result2 = get_mirna_mrna_pair(gff_db, region2, region1)

            if result2:
                miRNA, mRNA = result2
                features1_ids = ", ".join([feature.attributes["product"][0] for feature in miRNA])
                features2_ids = ", ".join([feature.attributes["gene"][0] for feature in mRNA])

                out.write(f"Region 1: {region2}, Region 2: {region1}; ")
                out.write(f"miRNA: {features1_ids}, mRNA: {features2_ids}\n")

                for miRNA_feature in miRNA:
                    for mRNA_feature in mRNA:
                        interaction_counts[(miRNA_feature.attributes["product"][0], mRNA_feature.attributes["gene"][0])] += 1

        out.write("\n")
        out.write("###############################################\n")

        found_interactions = 0
        not_found_interactions = 0
        not_found_miRNA = 0

        mirtar_interactions = load_mirtar_interactions(mirtar_file)

        for interaction, count in interaction_counts.items():
            out.write(f"{interaction}: {count} ")

            miRNA_id = interaction[0]

            if miRNA_id in mirtar_interactions:
                if interaction[1] in mirtar_interactions[miRNA_id]:
                    out.write("FOUND\n")
                    found_interactions += count
                else:
                    out.write("NOT_FOUND\n")
                    not_found_interactions += count
            else:
                # If nor found try to remove the additional identifier if exists (three -)
                if miRNA_id.count("-") > 2:
                    # Remove last part
                    miRNA_id = miRNA_id.rsplit("-", 1)[0]
                    if miRNA_id in mirtar_interactions:
                        if interaction[1] in mirtar_interactions[miRNA_id]:
                            out.write("FOUND_REDUCED\n")
                            found_interactions += count
                        else:
                            out.write("INTERACTION_NOT_FOUND\n")
                            not_found_interactions += count
                    else:
                        out.write("MIRNA_REDUCED_NOT_FOUND\n")
                        not_found_miRNA += count
                else:
                    out.write("MIRNA_NOT_FOUND\n")
                    not_found_interactions += count

        out.write(f"\nFound interactions: {found_interactions}\n")
        out.write(f"Not found interactions: {not_found_interactions}\n")
        out.write(f"Not found miRNA: {not_found_miRNA}\n")




if __name__ == "__main__":
    db_file = sys.argv[4]
    db = load_gff_db(sys.argv[1], db_file)
    bed_file = sys.argv[2]
    output_file = sys.argv[3]
    mirtar_file = sys.argv[5]
    process_bed_file(bed_file, db, output_file)
