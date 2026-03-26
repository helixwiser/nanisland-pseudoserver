#!/usr/bin/env python3


import argparse
import pandas as pd
import glob
import os
import math



def get_protein_name(uniprot_id, details_dict):
    return details_dict.get(uniprot_id, "Unknown")




def main():
    parser = argparse.ArgumentParser(description="Redefine representative protein by median abundance.")
    parser.add_argument("--input-pattern", required=True, help="Glob pattern for *_membercount.csv files")
    parser.add_argument("--detail-file", required=True, help="Path to merged_details_unique.tsv")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    parser.add_argument("--tissue-type", choices=["cancer", "normal", "both"], default="both",
                        help="Process only cancer, normal, or both tissue types")


    args = parser.parse_args()


    input_pattern = args.input_pattern
    detail_file = args.detail_file
    output_dir = args.output_dir
    tissue_filter = args.tissue_type


    os.makedirs(output_dir, exist_ok=True)


    # Load protein details
    details_df = pd.read_csv(detail_file, sep="\t", usecols=["uniprot_id", "protein_name"])
    details_dict = pd.Series(details_df.protein_name.values, index=details_df.uniprot_id.values).to_dict()


    # Read all files
    file_paths = glob.glob(input_pattern)
    if not file_paths:
        raise FileNotFoundError("No files matched pattern: " + input_pattern)


    data_list = []


    for file_path in file_paths:
        filename = os.path.basename(file_path)
        if "cancer" in filename:
            tissue_type = "cancer"
        elif "normal" in filename:
            tissue_type = "normal"
        else:
            continue
        if tissue_filter != "both" and tissue_type != tissue_filter:
            continue
        patient_id = filename.split("_")[0]
        df = pd.read_csv(file_path, sep=",")
        df = df.dropna(subset=["proteinfamily"])
        df["Patient_ID"] = patient_id
        df["tissue_type"] = tissue_type
        data_list.append(df)


    if not data_list:
        raise ValueError("No valid data loaded from input files.")


    combined_df = pd.concat(data_list, ignore_index=True)


    # Standardize proteinfamily: take first ID if comma-separated
    combined_df["proteinfamily"] = combined_df["proteinfamily"].astype(str).str.split(",").str[0].str.strip()


    # Group by Patient_ID and proteinfamily
    def get_most_abundant(row):
        row = row.sort_values(by="count", ascending=False)
        top_row = row.iloc[0]
        return pd.Series({
            "mostabundantMember": top_row["proteinid"],
            "mostabundantMember_count": top_row["count"],
            "totalProteinfamily_count": row["count"].sum()
        })


    summary_df = combined_df.groupby(["Patient_ID", "proteinfamily"]).apply(get_most_abundant).reset_index()
    summary_df = summary_df.merge(combined_df[["Patient_ID", "tissue_type"]].drop_duplicates(), on="Patient_ID")


    # Save output
    if tissue_filter in ["cancer", "both"]:
        cancer_df = summary_df[summary_df["tissue_type"] == "cancer"].copy()
        cancer_out = os.path.join(output_dir, "combined_cancerskin_skeleton_membercount.csv")
        cancer_df.to_csv(cancer_out, index=False, sep=",")


    if tissue_filter in ["normal", "both"]:
        normal_df = summary_df[summary_df["tissue_type"] == "normal"].copy()
        normal_out = os.path.join(output_dir, "combined_normalskin_skeleton_membercount.csv")
        normal_df.to_csv(normal_out, index=False, sep=",")


    # Process per proteinfamily
    grouped = combined_df.groupby("proteinfamily")


    # Original rep: the proteinfamily ID itself
    orig_rep = grouped["proteinfamily"].first()

#
#    # Redefined rep: protein with highest **median** count across samples
#    median_by_member = combined_df.groupby(["proteinfamily", "proteinid"])["count"].median()
#    idx = median_by_member.groupby("proteinfamily").idxmax()
#
#
#    proteinid_list = []
#    pfam_list = []
#    for pfam, key in idx.items():
#        if isinstance(key, tuple) and len(key) == 2:
#            proteinid = key[1]
#        else:
#            proteinid = "Unknown"
#        proteinid_list.append(proteinid)
#        pfam_list.append(pfam)
#
#
#    redefined_rep = pd.Series(proteinid_list, index=pfam_list)

   # === Redefined rep: use weighted score = median_count * log10(detection_count + 1) ===
    grouped = combined_df.groupby(["proteinfamily", "proteinid"])
    median_count = grouped["count"].median()
    detection_count = grouped.size()  # number of samples this member appears in




    scores = {}
    for (pfam, pid), med in median_count.items():
        det = detection_count[(pfam, pid)]
        #score = med * math.log10(1 + det)
        score = math.log10(1 + med) * det
        if pfam not in scores:
            scores[pfam] = []
        scores[pfam].append((pid, score))




    redefined_rep = {}
    for pfam, candidates in scores.items():
        best_pid = max(candidates, key=lambda x: x[1])[0]
        redefined_rep[pfam] = best_pid




    redefined_rep = pd.Series(redefined_rep)
    # Median total family count across patients
    total_by_sample = combined_df.groupby(["Patient_ID", "proteinfamily"])["count"].sum()
    median_family_count = total_by_sample.groupby("proteinfamily").median()


    # Build summary
    family_summary = pd.DataFrame({
        "proteinfamily": orig_rep.index,
        "orig_repProtein": orig_rep.values,
        "redefined_repProtein": redefined_rep.reindex(orig_rep.index).fillna("Unknown").values,
        "medianProteinfamily_count": median_family_count.reindex(orig_rep.index).fillna(0).values
    })


    # Add names
    family_summary["orig_repProteinname"] = family_summary["orig_repProtein"].apply(
        lambda x: get_protein_name(x, details_dict))
    family_summary["redefined_repProteinname"] = family_summary["redefined_repProtein"].apply(
        lambda x: get_protein_name(x, details_dict))


    # Final columns
    family_summary = family_summary[[
        "proteinfamily", "orig_repProtein", "orig_repProteinname",
        "redefined_repProtein", "redefined_repProteinname", "medianProteinfamily_count"
    ]]


    # Save family-level counts
    if tissue_filter in ["cancer", "both"]:
        cancer_family_count = cancer_df.groupby("proteinfamily")["totalProteinfamily_count"].sum().reset_index()
        cancer_family_count.to_csv(os.path.join(output_dir, "proteinfamily_cancerskin_count.csv"), index=False, sep=",")


    if tissue_filter in ["normal", "both"]:
        normal_family_count = normal_df.groupby("proteinfamily")["totalProteinfamily_count"].sum().reset_index()
        normal_family_count.to_csv(os.path.join(output_dir, "proteinfamily_normalskin_count.csv"), index=False, sep=",")


    # Save redefined summary
    family_summary.to_csv(os.path.join(output_dir, "proteinfamily_redefined_summary.csv"), index=False, sep=",")


    print("Processing complete. Outputs written to " + output_dir)




if __name__ == "__main__":
    main()