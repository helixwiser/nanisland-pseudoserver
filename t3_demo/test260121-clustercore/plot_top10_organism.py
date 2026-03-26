import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt
from collections import Counter




def main(summary_file, go_file, output_dir):
    # Load data
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Summary file not found: {summary_file}")
    if not os.path.exists(go_file):
        raise FileNotFoundError(f"GO file not found: {go_file}")


    summary_df = pd.read_csv(summary_file)
    go_df = pd.read_csv(go_file)


    # Extract redefined representative proteins
    proteins = summary_df["redefined_repProtein"].dropna().unique()


    # Filter GO table for these proteins
    filtered_go = go_df[go_df["uniprot_id"].isin(proteins)]


    if filtered_go.empty:
        print("No matching proteins found in GO file.")
        return


    # Drop rows where 'organism' is NaN or not a string
    filtered_go = filtered_go.dropna(subset=["organism"])
    filtered_go = filtered_go[filtered_go["organism"].apply(lambda x: isinstance(x, str))]


    # Count occurrences of each organism
    org_counter = Counter(filtered_go["organism"])
    top10 = org_counter.most_common(10)
    top10_orgs = [str(org) for org, _ in top10]
    top10_counts = [count for _, count in top10]


    # Reverse for descending order (largest on top)
    top10_orgs = top10_orgs[::-1]
    top10_counts = top10_counts[::-1]


    # Create figure with better size and DPI
    plt.figure(figsize=(12, 8), dpi=150)


    # Use a clean bar color and add edge
    bars = plt.barh(top10_orgs, top10_counts, color="#2E8B57", edgecolor="black", linewidth=0.8)


    # Add count labels at the end of each bar
    for i, (bar, count) in enumerate(zip(bars, top10_counts)):
        plt.text(count + max(top10_counts) * 0.01, i, str(count), va="center", ha="left", fontsize=10, fontweight="bold")


    # Styling
    plt.xlabel("Number of Skeleton Nodes", fontsize=12, fontweight="bold")
    plt.ylabel("Organism", fontsize=12, fontweight="bold")
    plt.title("Top 10 Organisms in Skeleton Nodes", fontsize=14, fontweight="bold", pad=20)


    # Adjust margins and layout
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)


    # Grid for better readability
    plt.grid(axis="x", linestyle="--", alpha=0.7, linewidth=0.7)


    # Save plot
    plot_path = os.path.join(output_dir, "top10_organism_skeleton_nodes.png")
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()
    print(f"Plot saved to {plot_path}")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate top 10 organism plot from skeleton nodes using full organism name.")
    parser.add_argument("--summary-file", required=True, help="Path to proteinfamily_redefined_summary.csv")
    parser.add_argument("--go-file", required=True, help="Path to mergedids_GOstats.csv")
    parser.add_argument("--output-dir", default=".", help="Output directory for plot")
    args = parser.parse_args()
    main(args.summary_file, args.go_file, args.output_dir)