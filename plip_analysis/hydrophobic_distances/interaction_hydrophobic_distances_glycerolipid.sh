#!/bin/bash

root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/Gylcerolipids"
output_file="hydrophobic_distances_glycerolipid.txt"
temp_dir="temp_hydro_dist"

mkdir -p "$temp_dir"
> "$output_file"

find "$root_dir" -type f -name "report.txt" | while read -r report_path; do
    awk -v temp_dir="$temp_dir" '
        BEGIN { in_hydro_section = 0; outfile = temp_dir "/hydrophobic.tmp" }

        /^\*\*Hydrophobic Interactions\*\*/ { in_hydro_section = 1; next }
        /^\*\*/ { in_hydro_section = 0 }

        in_hydro_section && /^\| RESNR/ { next }

        in_hydro_section && /^\|/ && !/RESNR/ {
            gsub(/\|/, "", $0);
            gsub(/[ \t]+/, " ", $0);
            split($0, fields, " ");
            dist = fields[7];
            if (dist ~ /^[0-9.]+$/) {
                print dist >> outfile
            }
        }
    ' "$report_path"
done

echo "Hydrophobic interaction distances across all files:" >> "$output_file"
sort -n "$temp_dir/hydrophobic.tmp" >> "$output_file"

rm -r "$temp_dir"
echo "Consolidated hydrophobic distances saved in: $output_file"
