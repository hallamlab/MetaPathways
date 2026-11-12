import csv
from collections import defaultdict

annotations = defaultdict(set)
with open("panc-annotations.csv", newline="") as stream:
    for row in csv.DictReader(stream):
        annotations[(row["sample_id"], row["orf_id"])].add(row["reference_db"])

with open("panto-pathway-genes.csv", newline="") as stream:
    genes = list(csv.DictReader(stream))

with open("panto-evidence.csv", "w", newline="") as stream:
    fields = ["sample_id", "entity_id", "pathway_id", "orf_id", "product",
              "databases_in_panc_export"]
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    for gene in genes:
        evidence = annotations[(gene["sample_id"], gene["orf_id"])]
        writer.writerow({**{key: gene[key] for key in fields[:-1]},
                         "databases_in_panc_export": ";".join(sorted(evidence))})
