#Reads the clean_pairs.txt file and identifies the minimum number of mergers needed to achieve the desired task. To simplify this, we can create a graph where each orthogroup is a node, and each pair of orthogroups represents an edge. Then, we can find the connected components in the graph, which will give us the groups of orthogroups that need to be merged together.
## Usage: python calculate_min_mergers.py > min_mergers.txt
def read_pairs(file_path):
    pairs = []
    with open(file_path, 'r') as pairs_file:
        lines = pairs_file.readlines()

    i = 0
    while i + 1 < len(lines):
        orthogroup1_line = lines[i].strip()
        orthogroup2_line = lines[i + 1].strip()

        orthogroup1 = orthogroup1_line.split(":")[1].strip()
        orthogroup2 = orthogroup2_line.split(":")[1].strip()

        pairs.append((orthogroup1, orthogroup2))
        i += 2

    return pairs


def create_graph(pairs):
    graph = {}
    for orthogroup1, orthogroup2 in pairs:
        if orthogroup1 not in graph:
            graph[orthogroup1] = set()
        if orthogroup2 not in graph:
            graph[orthogroup2] = set()

        graph[orthogroup1].add(orthogroup2)
        graph[orthogroup2].add(orthogroup1)

    return graph


def find_connected_components(graph):
    components = []
    visited = set()

    def dfs(node, component):
        visited.add(node)
        component.add(node)

        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, component)

    for node in graph:
        if node not in visited:
            component = set()
            dfs(node, component)
            components.append(component)

    return components


def main():
    cleaner_pairs_file = 'cleaner_pairs.txt'
    min_mergers_file = 'min_mergers.txt'

    # Read pairs from the clean_pairs.txt file
    pairs = read_pairs(clean_pairs_file)

    # Create a graph from the pairs
    graph = create_graph(pairs)

    # Find connected components in the graph
    connected_components = find_connected_components(graph)

    # Open output file for writing
    with open(min_mergers_file, "w") as out:
        # Write the minimum number of mergers needed
        out.write(f"Minimum number of mergers needed: {len(connected_components)}\n\n")

        # Write each connected component group
        for i, component in enumerate(connected_components):
            out.write(f"Group {i + 1}: {', '.join(component)}\n")


if __name__ == "__main__":
    main()
