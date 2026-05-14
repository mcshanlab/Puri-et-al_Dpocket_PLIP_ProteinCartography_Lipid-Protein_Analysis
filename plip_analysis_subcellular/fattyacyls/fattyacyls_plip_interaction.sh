#!/bin/bash

root_dir="/Volumes/GigiMurin/plip/Fatty_acyl"
output_file="fattyacyls_plip_interaction_counts.tsv"

echo -e "file\thydrophobic\thbond\tsaltbridge\tpistacking\tpication\thalogen\tmetal" > "$output_file"

find "$root_dir" -type f -name "report.txt" -print0 | while IFS= read -r -d '' report_path; do

awk -v file="$report_path" '
BEGIN {
    hydro=0; hb=0; sb=0; pi=0; pc=0; halo=0; metal=0;
    sec=""
}

# -----------------------------
# NORMALIZED SECTION DETECTION
# -----------------------------
/^\*\*/ {

    line = tolower($0)

    if (line ~ /hydrophobic/) sec="hydro"
    else if (line ~ /hydrogen bonds/) sec="hb"
    else if (line ~ /salt bridges/) sec="sb"
    else if (line ~ /pi[- ]?stack/) sec="pi"
    else if (line ~ /cation/) sec="pc"
    else if (line ~ /halogen/) sec="halo"
    else if (line ~ /metal/) sec="metal"
    else sec=""

    next
}

# -----------------------------
# DATA ROWS
# -----------------------------
/^\|/ {

    if (sec == "") next
    if ($0 ~ /RESNR/) next

    if (length($0) < 5) next

    if (sec=="hydro") hydro++
    else if (sec=="hb") hb++
    else if (sec=="sb") sb++
    else if (sec=="pi") pi++
    else if (sec=="pc") pc++
    else if (sec=="halo") halo++
    else if (sec=="metal") metal++
}

END {
    print file "\t" hydro "\t" hb "\t" sb "\t" pi "\t" pc "\t" halo "\t" metal
}
' "$report_path"

done >> "$output_file"

echo "Saved: $output_file"
