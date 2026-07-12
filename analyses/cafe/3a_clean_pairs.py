# This script takes the grep_similarity.txt summary, and eliminated duplicated comparisons. For example, if OG1 and OG2 are 80% similar, we do not also need OG2 vs. OG1 similarity to know that the two need to be merged. 

seen_pairs = set()

with open('grep_similarity.txt', 'r') as file:
    lines = file.readlines()

clean_lines = []
for i in range(0, len(lines), 3):
    orthogroup1_line = lines[i].strip()
    orthogroup2_line = lines[i + 1].strip()
    similarity_line = lines[i + 2].strip()

    orthogroup1 = orthogroup1_line.split(":")[1].strip()
    orthogroup2 = orthogroup2_line.split(":")[1].strip()

    pair = (orthogroup1, orthogroup2)
    reversed_pair = (orthogroup2, orthogroup1)

    if pair not in seen_pairs and reversed_pair not in seen_pairs:
        clean_lines.append(f"{orthogroup1_line}\n")
        clean_lines.append(f"{orthogroup2_line}\n")
        clean_lines.append(f"{similarity_line}\n")

        seen_pairs.add(pair)

with open('clean_pairs.txt', 'w') as file:
    file.writelines(clean_lines)
