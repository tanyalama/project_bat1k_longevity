import csv
import re
from collections import Counter

# Combined regex pattern for olfactory receptor genes
or_pattern = re.compile(r'^OR\d+(?:[A-Z]+\d*|\d+)$')

def extract_gene_symbols(row, fieldnames):
    """Extract unique gene symbols from all species columns."""
    genes = []
    for field in fieldnames[1:]:
        entries = row[field].split(', ')
        for entry in entries:
            parts = entry.split('.')
            if len(parts) > 2:
            # remove first and last parts, join the middle parts back with '.'
                gene = '.'.join(parts[1:-1])
                if gene:
                    genes.append(gene)
            else:
            # otherwise just take the first part
                gene = parts[0]
                if gene:
                    genes.append(gene)
    return genes

def filter_olfactory_rows(input_file, output_file, log_file):
    with open(input_file, 'r') as infile, \
         open(output_file, 'w', newline='') as outfile, \
         open(log_file, 'w', newline='') as logfile:

        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.DictWriter(outfile, delimiter='\t', fieldnames=reader.fieldnames)
        log_writer = csv.writer(logfile, delimiter='\t')
        
        writer.writeheader()
        log_writer.writerow(['Orthogroup', 'Olfactory_Percent', 'Unique_Genes', 'OR_Genes'])

        for row in reader:
            gene_symbols = extract_gene_symbols(row, reader.fieldnames)
            unique_genes = set(gene_symbols)
            or_genes = [g for g in unique_genes if or_pattern.match(g)]
            olfactory_percent = len(or_genes) / len(unique_genes) if unique_genes else 0

            if olfactory_percent > 0.5:
                log_writer.writerow([
                    row[reader.fieldnames[0]],
                    f"{olfactory_percent:.2f}",
                    "; ".join(sorted(unique_genes)),
                    "; ".join(sorted(or_genes))
                ])
            else:
                writer.writerow(row)

    print(f"Finished filtering olfactory orthogroups.")
    print(f"Filtered output → {output_file}")
    print(f"Log of removed rows → {log_file}")

if __name__ == "__main__":
    filter_olfactory_rows(
        input_file='filtered_merged_orthogroups.tsv',
        output_file='no_olfactory_filtered_merged_orthogroups.tsv',
        log_file='removed_olfactory_log.tsv'
    )
