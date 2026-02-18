#!/bin/bash

# Set working directory and output file
root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/Fatty_acyl"
output_file="AA_fattyacyl.txt"
temp_dir="temp_amino_counts"

# Create temp directory for intermediate files
mkdir -p "$temp_dir"
> "$output_file"

# Traverse all report.txt files and collect counts per interaction type
find "$root_dir" -type f -name "report.txt" | while read -r report_path; do
    awk -v temp_dir="$temp_dir" '
        /^\*\*/ {
            section_name = $0;
            gsub(/\*\*/, "", section_name);
            gsub(/^[ \t]+|[ \t]+$/, "", section_name);  # Trim whitespace
            current_file = temp_dir "/" section_name ".tmp";
            in_table = 0;
            next;
        }

        /^\| RESNR/ { in_table = 1; next }
        /^\*\*/ { in_table = 0 }

        in_table && /^\|/ && !/RESNR/ {
            gsub(/\|/, "", $0);
            split($0, fields, " ");
            res_type = fields[2];
            if (res_type ~ /^[A-Z]{3}$/) {
                print res_type >> current_file;
            }
        }
    ' "$report_path"
done

# Consolidate into one output file
for f in "$temp_dir"/*.tmp; do
    interaction_type=$(basename "$f" .tmp)
    echo "### $interaction_type" >> "$output_file"
    sort "$f" | uniq -c | sort -nr >> "$output_file"
    echo "" >> "$output_file"
done

# Clean up
rm -r "$temp_dir"

echo "Consolidated amino acid counts saved in: $output_file"

python3 plot_AAs_interaction.py
