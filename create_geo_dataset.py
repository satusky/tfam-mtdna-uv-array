import os
import argparse
import numpy as np
import pandas as pd


def main(ARGS):
    # Read in data
    all_stats = pd.read_csv(ARGS.data_file, index_col=[0,1], skipinitialspace=True, low_memory=False)

    # Init geo df
    geo_columns = ["sequence", "start_coord", "end_coord"]
    geo_df = all_stats[geo_columns].xs("ctrl_30", level=1).reset_index().copy()

    # Get all double-stranded z scores that have sequences (i.e. not controls)
    z_df = all_stats.dropna(subset="sequence").copy()
    z_df = z_df.loc[z_df["strands"] != 1, "norm_median_z"].unstack(level=1).reset_index()
    z_df.columns.name = None

    # Remove unused sequences
    z_df = z_df[["universal_top" not in str(n) for n in z_df["name"]]]
    z_df = z_df[["_o1_A" not in str(n) for n in z_df["name"]]]
    z_df = z_df[["_o2_A" not in str(n) for n in z_df["name"]]]
    z_df = z_df[["_o1_T" not in str(n) for n in z_df["name"]]]
    z_df = z_df[["_o2_T" not in str(n) for n in z_df["name"]]]
    z_df = z_df[["_o1_G" not in str(n) for n in z_df["name"]]]
    z_df = z_df[["_o2_G" not in str(n) for n in z_df["name"]]]

    # Merge z scores into geo df
    merged = geo_df.merge(z_df, on="name", how="right")
    merged["ID_REF"] = np.arange(1, len(merged) + 1)
    merged = merged.rename(
        columns= {
            "name": "ID",
            "sequence": "Sequence",
            "start_coord": "Mito_genome_start_coord",
            "end_coord": "Mito_genome_end_coord",
            "ctrl_30": "TFAM_30nM",
            "ctrl_300": "TFAM_300nM",
            "uv_30": "TFAM_30nM_UVC",
            "uv_300": "TFAM_300nM_UVC"
        }
    )
    merged = merged[[
        "ID_REF", "ID", "Sequence", "Mito_genome_start_coord", "Mito_genome_end_coord",
        "TFAM_30nM", "TFAM_300nM", "TFAM_30nM_UVC", "TFAM_300nM_UVC"
    ]]

    os.makedirs(ARGS.output_dir, exist_ok=True)
    merged.to_excel(os.path.join(ARGS.output_dir, ARGS.output_file), index=False)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("--data_file", type=str, default="data/tfam_array_all_stats.csv", help="Path to data csv file")
    PARSER.add_argument("--output_file", type=str, default="GEO_matrix.xlsx", help="Path to output excel file")
    PARSER.add_argument("--output_dir", type=str, default=".", help="Output directory")
    ARGS = PARSER.parse_args()

    main(ARGS)
