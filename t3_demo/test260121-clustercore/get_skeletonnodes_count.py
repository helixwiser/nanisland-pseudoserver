import pandas as pd
import os
import argparse


def parse_read_file_for_counts(uniprot_read_file):
    id_to_count = {}
    id_to_name = {}


    if not os.path.exists(uniprot_read_file):
        raise FileNotFoundError(f"Read file not found: {uniprot_read_file}")


    with open(uniprot_read_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or "sp|" not in line:
                continue
            start = line.find("sp|")
            if start == -1:
                continue
            protein_part = line[start:]
            parts = protein_part.split("|", 2)
            if len(parts) < 3:
                continue
            protein_id = parts[1].strip()
            protein_name_full = parts[2].split()[0].strip()
            id_to_count[protein_id] = id_to_count.get(protein_id, 0) + 1
            if protein_id not in id_to_name:
                id_to_name[protein_id] = protein_name_full


    print(f"Parsed {len(id_to_count)} unique proteins from read file.")
    return id_to_count, id_to_name


def main(skeleton_edge_file, cluster_file, uniprot_read_file, output_file):
    if not os.path.exists(skeleton_edge_file):
        raise FileNotFoundError(f"File not found: {skeleton_edge_file}")
    edges = pd.read_csv(skeleton_edge_file)


    family_rep_to_members = {}
    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            rep, member = parts
            if rep not in family_rep_to_members:
                family_rep_to_members[rep] = []
            family_rep_to_members[rep].append(member)


    id_to_count, id_to_name = parse_read_file_for_counts(uniprot_read_file)


    output_rows = []
    for _, row in edges.iterrows():
        source = row["source"]
        target = row["target"]
        source_rep = f"AF-{source}-F1-model_v6"
        target_rep = f"AF-{target}-F1-model_v6"
        source_members = family_rep_to_members.get(source_rep, [source_rep])
        target_members = family_rep_to_members.get(target_rep, [target_rep])


        for member in source_members:
            if not (member.startswith("AF-") and "-F1-model_v6" in member):
                continue
            protein_id = member.split("-")[1]
            count = id_to_count.get(protein_id, 0)
            name = id_to_name.get(protein_id, "Unknown")
            output_rows.append({
                "proteinid": protein_id,
                "count": str(count),
                "name": name,
                "proteinfamily": source
            })


        for member in target_members:
            if not (member.startswith("AF-") and "-F1-model_v6" in member):
                continue
            protein_id = member.split("-")[1]
            count = id_to_count.get(protein_id, 0)
            name = id_to_name.get(protein_id, "Unknown")
            output_rows.append({
                "proteinid": protein_id,
                "count": str(count),
                "name": name,
                "proteinfamily": target
            })


    result_df = pd.DataFrame(output_rows)
    agg_df = result_df.groupby("proteinid").agg(
        count=("count", "max"),
        name=("name", "first"),
        proteinfamily=("proteinfamily", lambda x: ", ".join(set(x)))
    ).reset_index()


    agg_df.columns = ["proteinid", "count", "name", "proteinfamily"]
    agg_df.to_csv(output_file, index=False)
    print(f"Output saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count skeleton edge member occurrences using UniProt read mapping.")
    parser.add_argument("--skeleton-edge", required=True, help="Path to skeleton_edges.csv")
    parser.add_argument("--cluster", required=True, help="Path to cluster TSV file")
    parser.add_argument("--uniprot-read", required=True, help="Path to read-to-UniProt file")
    parser.add_argument("--output", default="nodesmembercount.csv", help="Output CSV file")
    args = parser.parse_args()
    main(args.skeleton_edge, args.cluster, args.uniprot_read, args.output)
