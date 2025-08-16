#!/bin/bash
# fetch_fields.sh
# Retrieve the full list of available ENA API fields for "sequence" results
# This ensures that all queries in fetch.sh are based on valid, up-to-date fields.
# The output is stored in sequence_fields.tsv for documentation and reproducibility.

set -euo pipefail

curl -s -G 'https://www.ebi.ac.uk/ena/portal/api/returnFields' \
  --data-urlencode 'dataPortal=ena' \
  --data-urlencode 'result=sequence' \
  -o sequence_fields.tsv

echo "Field definitions saved to sequence_fields.tsv"
