import csv
from collections import defaultdict

molecular_functions = defaultdict(set)

def get_file(file_path, name, column):
    with open(file_path, mode='r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            gene = row['Gene']
            functions = row[column].split(' | ')
            for function in functions:
                if function:
                    molecular_functions[function.strip()].add(gene)

    with open(f'../data/gsea/{name}.csv', mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow([column, 'Genes'])
        
        for function, genes in molecular_functions.items():
            writer.writerow([function, ', '.join(genes)])

get_file('../data/genes_with_annotations.csv', 'Dictybase_Phenotype', 'Phenotypes')