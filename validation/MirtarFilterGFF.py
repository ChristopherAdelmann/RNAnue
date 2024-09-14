
def load_mirtar_ids(mirtar_file: str):
    mirtar_ids = []

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

            mirtar_ids.append(mirtar_id)
            mirtar_ids.append(mirna_id)

    return mirtar_ids


def filter_gff(gff_file_in: str, gff_file_out, mirtar_ids: list):

    with open(gff_file_in) as file, open(gff_file_out, "w") as out:
        for line in file:
            if line.startswith("#"):
                out.write(line)
            else:
                for mirtar_id in mirtar_ids:
                    if mirtar_id in line:
                        out.write(line)
                        break


mirtar_file = "mirtarbase.txt"
gff_file_in = "genes.gff"
gff_file_out = "genes_mirtar.gff"

mirtar_ids = load_mirtar_ids(mirtar_file)
filter_gff(gff_file_in, gff_file_out, mirtar_ids)
