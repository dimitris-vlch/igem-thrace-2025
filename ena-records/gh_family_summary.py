#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# Script: gh_family_summary.py
# Author: Δημήτρης
# Date: 2025-08-18
#
# Description:
#   This script parses ENA record TSV files for glycoside hydrolase (GH) families
#   (e.g. GH18, GH19, GH23, GH48, GH75, GH80). It extracts and summarizes:
#     - total number of records
#     - unique species count
#     - unique countries represented
#     - environmental samples
#     - distinct hosts
#
# Input:
#   TSV files located in the folder "ena-records/", named as ghXX.tsv
#
# Output:
#   Summary statistics printed to stdout for each GH family
#
# Example usage:
#   python3 gh_family_summary.py
#
# Notes:
#   Ensure TSV files follow the ENA metadata format with appropriate columns.
# ------------------------------------------------------------------------------

import pandas as pd
from pathlib import Path

# Folder containing TSV files
folder = Path(".")

# Loop through all GH TSV files
for file in folder.glob("gh*.tsv"):
    print(f"\n=== Family {file.stem.upper()} ===")
    
    # Load TSV file into pandas DataFrame
    df = pd.read_csv(file, sep="\t")
    
    # Total number of records
    print(f"Total records: {len(df)}")
    
    # Number of unique species
    if "scientific_name" in df.columns:
        unique_species = df["scientific_name"].nunique()
        print(f"Distinct species: {unique_species}")
    
    # Number of unique countries (if column exists)
    if "countrylocation" in df.columns:
        unique_countries = df["countrylocation"].dropna().nunique()
        print(f"Countries: {unique_countries}")
    
    # Check how many are marked as environmental samples
    if "environmental_sample" in df.columns:
        env_samples = df["environmental_sample"].astype(str).str.lower().eq("true").sum()
        print(f"Environmental samples: {env_samples}")
    
    # Count hosts (if present)
    if "host" in df.columns:
        hosts = df["host"].dropna().nunique()
        print(f"Distinct hosts: {hosts}")
