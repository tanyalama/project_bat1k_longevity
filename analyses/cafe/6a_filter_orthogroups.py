python!
import csv

def process_gene_name(gene):
    parts = gene.split('.')
    if len(parts) > 2:
        # remove first and last parts, join the middle parts back with '.'
        return '.'.join(parts[1:-1])
    else:
        # otherwise just take the first part
        return parts[0]

def count_genes_per_species(row, species_columns):
    species_gene_counts = {}
    all_genes = set()
    for species in species_columns:
        raw_genes = row[species].split(', ') if row[species].strip() else []
        processed_genes = [process_gene_name(g) for g in raw_genes]
        species_gene_counts[species] = len(raw_genes)            # total genes count (no uniqueness here)
        all_genes.update(processed_genes)                         # unique genes across species (processed)
    return species_gene_counts, len(all_genes)

def check_filters(species_gene_counts, total_unique_genes, max_genes_threshold=450, max_variation_threshold=100):
    max_genes = max(species_gene_counts.values())
    min_genes = min(species_gene_counts.values())
    variation = max_genes - min_genes
    reasons = []
    if max_genes > max_genes_threshold:
        reasons.append(f"species with >{max_genes_threshold} genes (max: {max_genes})")
    if variation > max_variation_threshold:
        reasons.append(f"variation > {max_variation_threshold} genes (variation: {variation})")
    if total_unique_genes > max_genes_threshold:
        reasons.append(f"total unique genes > {max_genes_threshold} (total: {total_unique_genes})")
    return reasons

def filter_orthogroups(input_tsv, output_tsv, log_file):
    with open(input_tsv, newline='') as infile, \
         open(output_tsv, 'w', newline='') as outfile, \
         open(log_file, 'w') as log:
        
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        
        species_columns = reader.fieldnames[1:]  # Assuming first column is orthogroup name
        
        removed_count = 0
        kept_count = 0
        
        for row in reader:
            species_gene_counts, total_unique_genes = count_genes_per_species(row, species_columns)
            
            reasons = check_filters(species_gene_counts, total_unique_genes)
            
            if reasons:
                removed_count += 1
                log.write(f"Removed {row[reader.fieldnames[0]]}: {', '.join(reasons)}\n")
            else:
                kept_count += 1
                writer.writerow(row)
        
        log.write(f"\nSummary:\nTotal orthogroups processed: {removed_count + kept_count}\n")
        log.write(f"Orthogroups removed: {removed_count}\n")
        log.write(f"Orthogroups kept: {kept_count}\n")

if __name__ == "__main__":
    input_tsv = 'cleaned_merged_orthogroups.tsv'  # Your input file
    output_tsv = 'filtered_merged_orthogroups.tsv'       # Output filtered orthogroups
    log_file = 'orthogroup_filter_log.txt'        # Log file
    
    filter_orthogroups(input_tsv, output_tsv, log_file)
