#!/bin/bash

# --------------------------------------
# User configuration
# --------------------------------------
root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/Fatty_acyl"
output_dir="fattyacyl_saltbridge_stats"
mkdir -p "$output_dir"

# Output files
out_carb="$output_dir/CARBOXYLATE_FILES.txt"
out_phos="$output_dir/PHOSPHATE_FILES.txt"

# Clear old content
> "$out_carb"
> "$out_phos"

# --------------------------------------
# Process each report.txt
# --------------------------------------
find "$root_dir" -type f -name "report.txt" | while read -r report_path; do

    carb_count=0
    phos_count=0

    awk '
        BEGIN { in_sb = 0 }

        # Enter Salt Bridges section
        /^\*\*Salt Bridges\*\*/ { in_sb = 1; next }

        # Leave Salt Bridges section if another section begins
        /^\*\*/ && !/^\*\*Salt Bridges\*\*/ { in_sb = 0 }

        # Skip header row
        in_sb && /^\| RESNR/ { next }

        # Parse rows
        in_sb && /^\|/ {
            gsub(/\|/, "", $0)
            gsub(/[ \t]+/, " ", $0)
            split($0, f, " ")

            lig_group = f[10]

            if (lig_group == "Carboxylate") carb_count++
            if (lig_group == "Phosphate")   phos_count++
        }

        END {
            printf("%d %d\n", carb_count, phos_count) > "/dev/stderr"
        }
    ' "$report_path" 2> counts_tmp.txt

    read carb_count phos_count < counts_tmp.txt
    rm -f counts_tmp.txt

    # Output file names + counts
    [[ "$carb_count" -gt 0 ]] && echo -e "${report_path}\t${carb_count}" >> "$out_carb"
    [[ "$phos_count" -gt 0 ]] && echo -e "${report_path}\t${phos_count}" >> "$out_phos"

done

echo "Carboxylate files written to: $out_carb"
echo "Phosphate files written to:   $out_phos"
