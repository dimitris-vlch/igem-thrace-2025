#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# gh_fasta_fetch.sh
# Retrieve accession-only TSVs from the ENA Portal API for selected chitinolytic 
# GH families, then fetch the corresponding FASTA sequences from the ENA Browser API.
#
# Description:
#   - For each GH family (GH18, GH19, GH23, GH48, GH75, GH80):
#       1) Fetch accession-only TSV (header: accession) from ENA Portal API
#       2) Read accessions and retrieve the corresponding FASTA sequences
#
# Output (in current directory):
#   gh18.tsv, gh19.tsv, gh23.tsv, gh48.tsv, gh75.tsv, gh80.tsv
#   gh18.fasta, gh19.fasta, gh23.fasta, gh48.fasta, gh75.fasta, gh80.fasta
#
# Author: iGEM Thrace 2025 - Dry Lab
# Repository: https://github.com/dimitris-vlch/igem-thrace-2025
# ------------------------------------------------------------------------------

set -Eeuo pipefail

# ENA endpoints
PORTAL="https://www.ebi.ac.uk/ena/portal/api/search"
BROWSER="https://www.ebi.ac.uk/ena/browser/api"
FAMILIES=(GH18 GH19 GH23 GH48 GH75 GH80)

for FAMILY in "${FAMILIES[@]}"; do
  stem="$(echo "$FAMILY" | tr '[:upper:]' '[:lower:]')"
  tsv="${stem}.tsv"
  fasta="${stem}.fasta"

  echo "[${FAMILY}] Fetching TSV -> ${tsv}"
  curl -fsS -G "$PORTAL" \
    --data-urlencode "result=sequence" \
    --data-urlencode "query=description=\"${FAMILY}\"" \
    --data-urlencode "fields=accession" \
    --data-urlencode "format=tsv" \
    --data-urlencode "limit=0" \
    -o "$tsv"

  # Extract accessions (skip header line)
  echo "[${FAMILY}] Preparing FASTA list"
  mapfile -t ACCS < <(tail -n +2 "$tsv" | sed 's/\r$//')
  (( ${#ACCS[@]} > 0 )) || { echo "[${FAMILY}] Empty TSV, skipping."; continue; }

  csv="$(printf "%s," "${ACCS[@]}")"; csv="${csv%,}"

  echo "[${FAMILY}] Downloading FASTA -> ${fasta}"
  url="${BROWSER}/fasta/${csv}?download=true"
  if ! curl -fsS "$url" -o "$fasta"; then
    echo "[${FAMILY}] FASTA download failed, retrying..."
    sleep 1
    curl -fsS "$url" -o "$fasta"
  fi

  # Remove empty lines (optional)
  sed -i '/^$/d' "$fasta" || true

  if grep -q '^>' "$fasta"; then
    echo "[${FAMILY}] Success: ${fasta}"
  else
    echo "[${FAMILY}] Warning: No FASTA headers found. IDs may not be sequence accessions."
  fi
done

echo "All done."
