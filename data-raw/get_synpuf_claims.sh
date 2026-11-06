#!/usr/bin/env bash
# Download the CMS 2008-2010 DE-SynPUF claims (20 samples, about 7.4 GB,
# discarded after filtering) and keep the claims of the Webster-Clark cohort
# in data-raw/synpuf_claims/*.csv. Needs data-raw/webster_clark/.
# Run from the repository root: bash data-raw/get_synpuf_claims.sh
set -eu
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT="$ROOT/data-raw/synpuf_claims"
TMP="${SYNPUF_TMP:-$OUT/tmp}"
mkdir -p "$OUT" "$TMP"
IDS="$OUT/cohort_ids.txt"
Rscript -e 'raw <- haven::read_sas("data-raw/webster_clark/Data/baseccwcohort.sas7bdat"); writeLines(raw$ID, commandArgs(TRUE)[1])' "$IDS"
CMS="https://www.cms.gov/research-statistics-data-and-systems/downloadable-public-use-files/synpufs/downloads"
DL="https://downloads.cms.gov/files"
for f in inpatient outpatient pde carrier; do : > "$OUT/$f.csv"; done
for s in $(seq 1 20); do
  for f in inpatient outpatient pde carrierA carrierB; do
    case $f in
      inpatient)  u="$CMS/de1_0_2008_to_2010_inpatient_claims_sample_$s.zip" ;;
      outpatient) u="$CMS/de1_0_2008_to_2010_outpatient_claims_sample_$s.zip" ;;
      pde)        u="$DL/DE1_0_2008_to_2010_Prescription_Drug_Events_Sample_$s.zip" ;;
      carrierA)   u="$DL/DE1_0_2008_to_2010_Carrier_Claims_Sample_${s}A.zip" ;;
      carrierB)   u="$DL/DE1_0_2008_to_2010_Carrier_Claims_Sample_${s}B.zip" ;;
    esac
    out="$OUT/${f%[AB]}.csv"
    zip="$TMP/${f}_$s.zip"
    curl -sSL --retry 3 --max-time 1200 -o "$zip" "$u"
    if [ ! -s "$out" ]; then unzip -p "$zip" | head -1 | tr -d "\r" > "$out" || true; fi
    n=$(unzip -p "$zip" | LC_ALL=C awk -F, 'NR == FNR { ids[$1] = 1; next } ($1 in ids)' "$IDS" - | tr -d "\r" | tee -a "$out" | wc -l)
    echo "sample $s $f: $n rows $(date +%T)"
    rm -f "$zip"
  done
done
echo "DONE $(date)"
