#!/bin/bash
# ------------------------------------------------------------------------------
# fetch.sh
# Automated metadata retrieval from ENA Portal API for chitinolytic GH families
#
# Description:
#   This script queries the ENA Portal API for coding sequences (CDS) associated 
#   with six glycoside hydrolase (GH) families (GH18, GH19, GH23, GH48, GH75, GH80).
#   For each GH family, a tab-delimited file is generated containing sequence 
#   accessions and associated metadata fields of ecological and phylogenetic interest.
#
# Metadata fields included:
#   accession, description, scientific_name, tax_id, tax_lineage, country, location, 
#   altitude, isolation_source, environmental_sample, host, host_tax_id, dev_stage, 
#   sex, tissue_type, disease, collected_by, identified_by, collection_date, 
#   collection_date_start, collection_date_end, isolate, strain, sub_species, 
#   sub_strain, cultivar, ecotype, variety, specimen_voucher, culture_collection, 
#   bio_material, study_accession, sample_accession, secondary_sample_accession
#
# Output:
#   One TSV file per GH family (e.g., gh18.tsv, gh19.tsv, … gh80.tsv).
#
# Author: iGEM Thrace 2025 - Dry Lab
# Repository: https://github.com/dimitris-vlch/igem-thrace-2025
# ------------------------------------------------------------------------------

# ENA API base URL
BASE_URL="https://www.ebi.ac.uk/ena/portal/api/search"

# Common fields to retrieve (must exist in ENA documentation)
FIELDS="accession,description,scientific_name,tax_id,tax_lineage,country,location,altitude,isolation_source,environmental_sample,host,host_tax_id,dev_stage,sex,tissue_type,collected_by,identified_by,collection_date,collection_date_start,collection_date_end,isolate,strain,sub_species,sub_strain,cultivar,ecotype,variety,specimen_voucher,culture_collection,bio_material,study_accession,sample_accession,secondary_sample_accession"

# GH families of interest
FAMILIES=("GH18" "GH19" "GH23" "GH48" "GH75" "GH80")

# Loop through each family and fetch data
for FAMILY in "${FAMILIES[@]}"; do
    echo "Fetching records for $FAMILY..."
    curl -s -G "$BASE_URL" \
        --data-urlencode "result=sequence" \
        --data-urlencode "query=description=\"$FAMILY\"" \
        --data-urlencode "fields=$FIELDS" \
        --data-urlencode "format=tsv" \
        -o "${FAMILY,,}.tsv"
    echo "Saved to ${FAMILY,,}.tsv"
done

echo "All GH family datasets retrieved successfully."

