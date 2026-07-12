import csv
import re

# === Config ===
orthogroups_file = 'Orthogroups.tsv'
min_mergers_file = 'min_mergers.txt'
merged_output_file = 'merged_output.tsv'
non_merged_output_file = 'non_merged_output.tsv'
combined_output_file = 'cleaned_merged_orthogroups.tsv'
summary_log_file = 'summary_mergers.txt'

# === Helper: natural sort key ===
def natural_key(text):
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', text)]

# === Read header ===
with open(orthogroups_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    header = reader.fieldnames

# === Parse mergers ===
merged_rows = {}
merged_orthogroups_set = set()
summary_log = []

group_counter = 1

with open(min_mergers_file, 'r') as f:
    for line in f:
        if line.startswith("Group"):
            _, group_data = line.strip().split(":")
            orthogroups = [og.strip() for og in group_data.split(",")]
            group_name = f"og_{group_counter}"
            group_counter += 1
            merged_row = {col: "" for col in header}
            merged_row[header[0]] = group_name

            for og in orthogroups:
                merged_orthogroups_set.add(og)
                with open(orthogroups_file, 'r') as cf:
                    reader = csv.DictReader(cf, delimiter='\t')
                    for row in reader:
                        if row[header[0]] == og:
                            for col in header[1:]:
                                values = row[col].split(',') if row[col] else []
                                existing = merged_row[col].split(',') if merged_row[col] else []
                                merged = existing + values
                                # Clean and deduplicate
                                cleaned = [v.strip() for v in merged if v.strip()]
                                deduped = sorted(set(cleaned), key=str)
                                merged_row[col] = ', '.join(deduped)
                            break
            merged_rows[group_name] = merged_row
            summary_log.append(f"{group_name}: {', '.join(orthogroups)}")

# === Parse original Orthogroups.tsv and find non-merged ===
non_merged_rows = {}
with open(orthogroups_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        og_name = row[header[0]]
        if og_name not in merged_orthogroups_set:
            # Clean each cell
            for col in header[1:]:
                values = row[col].split(',') if row[col] else []
                cleaned = [v.strip() for v in values if v.strip()]
                row[col] = ', '.join(sorted(set(cleaned), key=str))
            non_merged_rows[og_name] = row

# === Write merged_output.tsv ===
with open(merged_output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=header)
    writer.writeheader()
    for row in sorted(merged_rows.values(), key=lambda x: natural_key(x[header[0]])):
        writer.writerow(row)

# === Write non_merged_output.tsv ===
with open(non_merged_output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=header)
    writer.writeheader()
    for row in sorted(non_merged_rows.values(), key=lambda x: natural_key(x[header[0]])):
        writer.writerow(row)

# === Write cleaned_merged_orthogroups.tsv ===
with open(combined_output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=header)
    writer.writeheader()
    combined = list(merged_rows.values()) + list(non_merged_rows.values())
    for row in sorted(combined, key=lambda x: natural_key(x[header[0]])):
        writer.writerow(row)

# === Write summary_mergers.txt ===
with open(summary_log_file, 'w') as f:
    f.write("Merged Groups Summary:\n")
    for entry in summary_log:
        f.write(f"{entry}\n")

print("All outputs generated:")
print(f" - {merged_output_file}")
print(f" - {non_merged_output_file}")
print(f" - {combined_output_file}")
print(f" - {summary_log_file}")
