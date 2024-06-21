import csv

FILE = "Dictybase-Phenotypes"
def parse_line(line):
    parts = line.strip().split('\t')
    go_id = parts[0]
    description= parts[1].rsplit(',', 1)
    description = description[0].split(',', 3)[3]
    description = description[:-4]
    genes = parts[2:]
    return go_id, description, genes

with open(f'../data/genesets/{FILE}.gmt', 'r') as file:
    lines = file.readlines()

with open(f'../data/genesets/{FILE}_mod.gmt', 'w') as gmt_file:
    for line in lines:
        go_id, description, genes = parse_line(line)
        gmt_file.write(f"{go_id}\t\"{description}\"\t" + "\t".join(genes) + "\n")

print("GMT file has been created.")
