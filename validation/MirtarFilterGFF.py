
def load_mirtar_ids(mirtar_file: str):
    mirtar_ids = set()

    with open(mirtar_file) as file:
        file.readline()
        for line in file:
            (
                mirtar_id,
                mirna_id,
                mirna_species,
                target_gene_id,
                target_gene_entrez_id,
                target_gene_species,
                experiment_type,
                support_type,
                reference,
            ) = line.strip().split("\t")

            mirtar_ids.add(mirtar_id)
            mirtar_ids.add(mirna_id)

    return mirtar_ids

def get_gff_attributes(attributes: str) -> dict:
    attribute_dict = {}
    for attribute in attributes.split(";"):
        key, value = attribute.split("=")
        attribute_dict[key] = value
    return attribute_dict

def filter_gff(gff_file_in: str, gff_file_out, mirtar_ids: set):

    with open(gff_file_in) as file, open(gff_file_out, "w") as out:
        processed = 0
        included = 0
        for line in file:
            processed += 1
            if line.startswith("#"):
                out.write(line)
            else:
                (
                    chromosome,
                    source,
                    feature_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attributes,
                ) = line.strip().split("\t")

                if feature_type == "gene" or feature_type == "miRNA":
                    attributes = get_gff_attributes(attributes)

                    if "ID" in attributes:
                        gene_id = attributes["ID"]
                        if gene_id in mirtar_ids:
                            out.write(line)
                            included += 1
                            continue
                    if "product" in attributes:
                        gene_id = attributes["product"]
                        if gene_id in mirtar_ids:
                            out.write(line)
                            included += 1
                            continue


mirtar_file = "mirtarbase.txt"
gff_file_in = "genes.gff"
gff_file_out = "genes_mirtar.gff"

mirtar_ids = load_mirtar_ids(mirtar_file)
filter_gff(gff_file_in, gff_file_out, mirtar_ids)
