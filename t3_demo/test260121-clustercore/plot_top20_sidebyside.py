import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt
from collections import Counter




def load_protein_list(summary_file):
    """Load redefined representative proteins from summary file."""
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"File not found: {summary_file}")
    df = pd.read_csv(summary_file)
    return df["redefined_repProtein"].dropna().unique().tolist()




def get_organism_counts(proteins, go_file):
    """Map proteins to organism and count occurrences."""
    if not os.path.exists(go_file):
        raise FileNotFoundError(f"GO file not found: {go_file}")
    go_df = pd.read_csv(go_file)


    # Filter and clean
    filtered = go_df[go_df["uniprot_id"].isin(proteins)]
    filtered = filtered.dropna(subset=["organism"])
    filtered = filtered[filtered["organism"].apply(lambda x: isinstance(x, str))]


    return Counter(filtered["organism"])




def main(ca_summary, no_summary, go_file, output_dir):
    # Load protein lists
    ca_proteins = load_protein_list(ca_summary)
    no_proteins = load_protein_list(no_summary)


    # Get organism counts
    ca_orgs = get_organism_counts(ca_proteins, go_file)
    no_orgs = get_organism_counts(no_proteins, go_file)


    # Combine all organisms and get top 20
    all_orgs = set(ca_orgs.keys()) | set(no_orgs.keys())
    total_counts = {org: ca_orgs.get(org, 0) + no_orgs.get(org, 0) for org in all_orgs}
    top20_orgs = sorted(total_counts, key=lambda x: total_counts[x], reverse=True)[:20]


    # Prepare data
    ca_counts = [ca_orgs.get(org, 0) for org in top20_orgs]
    no_counts = [no_orgs.get(org, 0) for org in top20_orgs]


    # Reverse for descending order (largest on top)
    top20_orgs = top20_orgs[::-1]
    ca_counts = ca_counts[::-1]
    no_counts = no_counts[::-1]


    # Plot
    plt.figure(figsize=(14, 10), dpi=150)
    y_pos = range(len(top20_orgs))
    bar_height = 0.4


    plt.barh(y_pos, ca_counts, height=bar_height, color="#2E8B57", edgecolor="black", label="CA Skeleton", linewidth=0.8)
    plt.barh([y + bar_height for y in y_pos], no_counts, height=bar_height, color="#CD5C5C", edgecolor="black", label="NO Skeleton", linewidth=0.8)


    # Labels
    for i, (ca_val, no_val) in enumerate(zip(ca_counts, no_counts)):
        if ca_val > 0:
            plt.text(ca_val + max(ca_counts) * 0.01, i, str(ca_val), va="center", ha="left", fontsize=9, fontweight="bold")
        if no_val > 0:
            plt.text(no_val + max(no_counts) * 0.01, i + bar_height, str(no_val), va="center", ha="left", fontsize=9, fontweight="bold")


    # Styling
    plt.yticks([y + bar_height/2 for y in y_pos], [str(org) for org in top20_orgs], fontsize=10)
    plt.xlabel("Number of Skeleton Nodes", fontsize=12, fontweight="bold")
    plt.ylabel("Organism", fontsize=12, fontweight="bold")
    plt.title("Top 20 Organisms in CA vs NO Skeleton Nodes", fontsize=14, fontweight="bold", pad=20)
    plt.legend(loc="lower right", fontsize=11)
    plt.grid(axis="x", linestyle="--", alpha=0.6, linewidth=0.6)
    plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)


    # Save
    output_path = os.path.join(output_dir, "top20_organism_ca_vs_no_skeleton.png")
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()
    print(f"Combined plot saved to {output_path}")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare top 20 organisms in CA and NO skeleton nodes.")
    parser.add_argument("--ca-summary", required=True, help="Path to CA proteinfamily_redefined_summary.csv")
    parser.add_argument("--no-summary", required=True, help="Path to NO proteinfamily_redefined_summary.csv")
    parser.add_argument("--go-file", required=True, help="Path to mergedids_GOstats.csv")
    parser.add_argument("--output-dir", default=".", help="Output directory for plot")
    args = parser.parse_args()
    main(args.ca_summary, args.no_summary, args.go_file, args.output_dir)