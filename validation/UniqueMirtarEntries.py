import sys
from collections import defaultdict

def count_unique_mirtar_entries(mirtar_file: str):
    mirtar_ids = defaultdict(int)
    with open(mirtar_file) as file:
        for line in file:
            id = line.strip().split("\t")[20].split(";")[0].split("=")[1]
            mirtar_ids[id] += 1


    threshold = 5

    bigger_than_threshold = 0
    smaller_than_threshold = 0

    for id, count in mirtar_ids.items():
        if count > threshold:
            print(f"{id}: {count}")
            bigger_than_threshold += 1
        else:
            smaller_than_threshold += 1

    print(f"Total: {len(mirtar_ids)}")
    print(f"Entries bigger than {threshold}: {bigger_than_threshold}")
    print(f"Entries smaller than {threshold}: {smaller_than_threshold}")



if __name__ == "__main__":
    mirtar_file = sys.argv[1]
    count_unique_mirtar_entries(mirtar_file)
