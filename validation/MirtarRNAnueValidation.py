from dataclasses import dataclass
from enum import Enum
import sys
from collections import defaultdict
import re


@dataclass
class MirtarInteraction:
    mirtar_id: str
    mirna_id: str
    mirna_species: str
    target_gene_id: str
    target_gene_entrez_id: int
    target_gene_species: str
    experiment_type: str
    support_type: str
    reference: str


@dataclass
class RNAnueInteraction:
    cluster_id: str
    miRNA_feature_id: str
    target_gene_feature_id: str
    target_gene_entrez_id: int


@dataclass
class MiRNAFeature:
    id: str
    chromosome: str
    strand: str
    start: int
    end: int
    product_name: str


@dataclass
class ExonFeature:
    id: str
    chromosome: str
    strand: str
    start: int
    end: int
    gene_id: int


def interactions_by_miRNA_from_database_file(
    file_path: str,
) -> dict[str, list[MirtarInteraction]]:
    interactions_by_mirna_id = defaultdict(list)

    with open(file_path) as file:
        # Skip the header
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

            interactions_by_mirna_id[mirna_id].append(
                MirtarInteraction(
                    mirtar_id,
                    mirna_id,
                    mirna_species,
                    target_gene_id,
                    int(target_gene_entrez_id),
                    target_gene_species,
                    experiment_type,
                    support_type,
                    reference,
                )
            )
    return interactions_by_mirna_id


def get_gff_attributes(attributes: str) -> dict:
    attribute_dict = {}
    for attribute in attributes.split(";"):
        key, value = attribute.split("=")
        attribute_dict[key] = value
    return attribute_dict


@dataclass
class GFFFeatures:
    exon_features: dict[str, ExonFeature]
    miRNA_features: dict[str, MiRNAFeature]


def read_feature_gff_file(
    file_path: str,
) -> GFFFeatures:
    print("Reading miRNA features from GFF file: ", file_path)

    mirna_features_by_id = dict()
    exon_features_by_id = dict()

    with open(file_path) as file:
        for line in file:
            if line.startswith("#"):
                continue
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
                if feature_type == "miRNA":
                    attribute_dict = get_gff_attributes(attributes)
                    name = attribute_dict.get("product", None)

                    if not name:
                        print(
                            "No product name found for miRNA feature: ",
                            attribute_dict["ID"],
                            " line: ",
                            line.strip(),
                        )
                        continue

                    print("ID: ", attribute_dict["ID"])

                    mirna_features_by_id[attribute_dict["ID"]] = MiRNAFeature(
                        attribute_dict["ID"],
                        chromosome,
                        strand,
                        int(start),
                        int(end),
                        name,
                    )
                elif feature_type == "exon":
                    attribute_dict = get_gff_attributes(attributes)
                    dbxref = attribute_dict.get("Dbxref", None)

                    if dbxref:
                        refs = dbxref.split(",")
                        ref_dict = {}
                        for ref in refs:
                            values = ref.split(":")

                            if len(values) != 2:
                                continue

                            ref_dict[values[0]] = values[1]

                        gene_id = ref_dict.get("GeneID", None)

                        if gene_id:
                            exon_features_by_id[attribute_dict["ID"]] = ExonFeature(
                                attribute_dict["ID"],
                                chromosome,
                                strand,
                                int(start),
                                int(end),
                                gene_id,
                            )

                            continue

                    print("No gene ID found for exon feature: ", attribute_dict["ID"])

    print(
        "Found ",
        len(mirna_features_by_id),
        " miRNA features and , ",
        len(exon_features_by_id),
        " exon features in GFF file",
    )
    return GFFFeatures(
        miRNA_features=mirna_features_by_id, exon_features=exon_features_by_id
    )


def find_feature_by_id(
    id: str, gff_features: GFFFeatures
) -> MiRNAFeature | ExonFeature | None:

    if id in gff_features.miRNA_features:
        return gff_features.miRNA_features[id]

    if id in gff_features.exon_features:
        return gff_features.exon_features[id]

    return None


class feature_type(Enum):
    protein_coding = 1
    unknown = 2
    miRNA = 3


def get_feature_type(
    feature_id: str,
    gff_features: GFFFeatures,
) -> feature_type:

    feature = find_feature_by_id(feature_id, gff_features)

    if feature:
        if isinstance(feature, MiRNAFeature):
            return feature_type.miRNA

        if isinstance(feature, ExonFeature):
            return feature_type.protein_coding

        return feature_type.unknown

    return feature_type.unknown


def load_miRNA_interactions_from_rnanue_file(
    file_path: str, gff_features: GFFFeatures
) -> list[RNAnueInteraction]:

    print("Filtering miRNA interactions from RNAnue file: ", file_path)

    protein_codings = 0
    miRNAs = 0
    unknowns = 0

    feature_not_found = 0
    supplementary_feature_not_found = 0

    interactions = []
    with open(file_path) as file:
        for line in file:
            (
                cluster_id,
                first_feature_id,
                first_segment_chromosome,
                first_segment_strand,
                first_segment_start,
                first_segment_end,
                second_feature_id,
                second_segment_chromosome,
                second_segment_strand,
                second_segment_start,
                second_segment_end,
                number_of_reads,
                gcs,
                ghs,
                p_value,
                p_value_adjusted,
            ) = line.strip().split("\t")

            first_feature = find_feature_by_id(first_feature_id, gff_features)
            second_feature = find_feature_by_id(second_feature_id, gff_features)

            pattern = re.compile(
                r"\b[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\b"
            )

            if not first_feature:
                if not pattern.match(first_feature_id):
                    feature_not_found += 1
                    print(
                        "First feature not found: ",
                        first_feature_id,
                        " Interaction id: ",
                        cluster_id,
                    )
                    continue

                supplementary_feature_not_found += 1
                continue

            if not second_feature:
                if not pattern.match(second_feature_id):
                    feature_not_found += 1
                    print(
                        "Second feature not found: ",
                        second_feature_id,
                        " Interaction id: ",
                        cluster_id,
                    )
                    continue

                supplementary_feature_not_found += 1
                continue

            miRNA_id: str | None = None

            target_gene_entrez_id: int | None = None
            target_gene_id: str | None = None

            if isinstance(first_feature, MiRNAFeature):
                miRNAs += 1
                miRNA_id = first_feature.product_name
            elif isinstance(first_feature, ExonFeature):
                protein_codings += 1
                target_gene_entrez_id = first_feature.gene_id
                target_gene_id = first_feature.id
            else:
                unknowns += 1

            if isinstance(second_feature, MiRNAFeature):
                miRNAs += 1
                miRNA_id = second_feature.product_name
            elif isinstance(second_feature, ExonFeature):
                protein_codings += 1
                target_gene_entrez_id = second_feature.gene_id
                target_gene_id = second_feature.id
            else:
                unknowns += 1

            if miRNA_id and target_gene_entrez_id and target_gene_id:
                interactions.append(
                    RNAnueInteraction(
                        cluster_id, miRNA_id, target_gene_id, target_gene_entrez_id
                    )
                )

    print(
        f"Found {miRNAs} miRNAs, {protein_codings} protein codings and {unknowns} unknown features in RNAnue file"
    )
    print(
        "Found ",
        len(interactions),
        " miRNA interactions in RNAnue file, missed ",
        supplementary_feature_not_found,
        " supplementary features and ",
        feature_not_found,
        " other features.",
    )
    return interactions


def find_miRNA_interaction_by_product_name(
    product_name: str, miRNA_features: list[MiRNAFeature]
) -> bool:
    for feature in miRNA_features:
        if feature.product_name == product_name:
            return True

    return False


def validate_gff_features_against_mirtar_database(
    mirtar_interactions: dict[str, list[MirtarInteraction]],
    miRNA_features: dict[str, MiRNAFeature],
):
    found = 0
    not_found = 0

    for interactions in mirtar_interactions.values():
        for interaction in interactions:
            if not find_miRNA_interaction_by_product_name(
                interaction.mirna_id, list(miRNA_features.values())
            ):
                not_found += 1
            else:
                found += 1

    print(
        "Found ",
        found,
        " miRNA from mirtar interactions in miRNA features from gff file",
    )
    print(
        "Not found ",
        not_found,
        " miRNA from mirtar interactions in miRNA features from gff file",
    )


def search_mirtar_interactions_by_rnanue_interaction(
    rnanue_interaction: RNAnueInteraction,
    mirtar_interactions: dict[str, list[MirtarInteraction]],
) -> list[MirtarInteraction] | None:

    if rnanue_interaction.miRNA_feature_id in mirtar_interactions:
        return mirtar_interactions[rnanue_interaction.miRNA_feature_id]

    return None


def validate_rnanue_interactions_against_mirtar_database(
    rnanue_interactions: list[RNAnueInteraction],
    mirtar_interactions: dict[str, list[MirtarInteraction]],
):
    assert isinstance(rnanue_interactions[0], RNAnueInteraction)

    found = 0
    not_found = 0

    for interaction in rnanue_interactions:
        database_results = search_mirtar_interactions_by_rnanue_interaction(
            interaction, mirtar_interactions
        )

        if not database_results:
            not_found += 1
        else:
            found += 1


def main():
    # Load command line arguments #1 mirtar file, #2 miRNA gff file, #3 rnanue file
    mirtar_file_path = sys.argv[1]
    gff_file_path = sys.argv[2]
    rnanue_file_path = sys.argv[3]

    gff_features = read_feature_gff_file(gff_file_path)

    mirtar_interactions = interactions_by_miRNA_from_database_file(mirtar_file_path)

    validate_gff_features_against_mirtar_database(
        mirtar_interactions, gff_features.miRNA_features
    )

    rnanue_mirna_interactions = load_miRNA_interactions_from_rnanue_file(
        rnanue_file_path, gff_features
    )

    if len(rnanue_mirna_interactions) == 0:
        print("No miRNA interactions found in RNAnue file")
        return

    validate_rnanue_interactions_against_mirtar_database(
        rnanue_mirna_interactions, mirtar_interactions
    )


if __name__ == "__main__":
    main()
