import sys
import os
from Bio import SeqIO

def load_mapping(tsv_path):
    mapping = {}
    with open(tsv_path) as f:
        header = next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            enst_id, gene_name = parts[0], parts[1]
            base_id = enst_id.split('.')[0]
            mapping[base_id] = gene_name
    return mapping

def strip_version(enst_id):
    return enst_id.split('.')[0]

def get_version(enst_id):
    return enst_id.split('.')[1] if '.' in enst_id else ""

def annotate_fasta(fasta_file, mapping):
    out_file = os.path.splitext(fasta_file)[0] + "_annotated.faa"
    print(f"Annotating {fasta_file} → {out_file}")
    with open(out_file, "w") as out_f:
        for rec in SeqIO.parse(fasta_file, "fasta"):
            base_id = strip_version(rec.id)
            version = get_version(rec.id)
            gene_symbol = mapping.get(base_id)
            if gene_symbol:
                rec.id = f"{base_id}.{gene_symbol}.{version}" if version else f"{base_id}.{gene_symbol}"
            else:
                rec.id = rec.id
            rec.description = ""
            rec.name = ""
            SeqIO.write(rec, out_f, "fasta")
    print(f"Done writing {out_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python get_faas_gene_names.py transcript_to_symbol.tsv list_of_faa_files.txt")
        sys.exit(1)

    mapping_file = sys.argv[1]
    list_file = sys.argv[2]

    mapping = load_mapping(mapping_file)
    print(f"Loaded {len(mapping)} transcript→gene mappings")

    with open(list_file) as f:
        fasta_files = [line.strip() for line in f if line.strip()]

    for fasta in fasta_files:
        if not os.path.isfile(fasta):
            print(f"Warning: file {fasta} not found, skipping")
            continue
        annotate_fasta(fasta, mapping)

if __name__ == "__main__":
    main()
