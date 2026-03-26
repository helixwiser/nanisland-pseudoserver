#!/usr/bin/env python3


import pandas as pd
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Annotate hub nodes with redefined protein family info.")
    parser.add_argument("--hub-file", required=True, help="Path to ppi_hub_nodes_annotated.csv")
    parser.add_argument("--summary-file", required=True, help="Path to proteinfamily_redefined_summary.csv")
    parser.add_argument("--details-file", required=True, help="Path to merged_details_unique.tsv")
    parser.add_argument("--output-file", required=True, help="Output path for annotated hub nodes")


    args = parser.parse_args()


    try:
        # Load data
        ppi_hub = pd.read_csv(args.hub_file)
        proteinfamily_summary = pd.read_csv(args.summary_file)
        merged_details = pd.read_csv(args.details_file, sep="\t", header=0)


        # Rename columns for clarity (ensure consistency)
        merged_details.rename(columns={
            "uniprot_id": "uniprot_id",
            "gene_name": "gene_name",
            "protein_name": "protein_name",
            "Organism": "Organism"
        }, inplace=True)


        ppi_hub.rename(columns={"Node": "Node"}, inplace=True)
        proteinfamily_summary.rename(columns={
            "proteinfamily": "proteinfamily",
            "orig_repProtein": "orig_repProtein",
            "redefined_repProtein": "redefined_repProtein",
            "redefined_repProteinname": "redefined_repProteinname",
            "medianProteinfamily_count": "medianProteinfamily_count"
        }, inplace=True)


        # Step 1: Merge hub nodes with protein family summary
        merged = pd.merge(ppi_hub, proteinfamily_summary, left_on="Node", right_on="orig_repProtein", how="left")


        # Step 2: Prepare details table for merging using redefined_repProtein
        details_for_redefined = merged_details[["uniprot_id", "gene_name", "protein_name", "Organism"]].copy()
        details_for_redefined.rename(columns={"uniprot_id": "redefined_repProtein"}, inplace=True)


        # Merge to get gene name, protein name, and organism for redefined protein
        merged = pd.merge(merged, details_for_redefined, on="redefined_repProtein", how="left")


        # Rename final columns
        merged.rename(columns={
            "gene_name": "redefined_repGeneName",
            "protein_name": "redefined_repProteinName"
        }, inplace=True)


        # Select final columns
        output_columns = [
            "Node",
            "orig_repProtein",
            "redefined_repProtein",
            "medianProteinfamily_count",
            "redefined_repGeneName",
            "redefined_repProteinName",
            "Organism",
            "Degree",
            "Degree_Centrality",
            "Betweenness",
            "Closeness",
            "PageRank"
        ]


        # Fill NaN with empty string (as per requirement)
        result = merged[output_columns].fillna("")


        # Save result
        result.to_csv(args.output_file, index=False)


        print(f"Success: Annotated hub nodes saved to {args.output_file}")
        return 0


    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())