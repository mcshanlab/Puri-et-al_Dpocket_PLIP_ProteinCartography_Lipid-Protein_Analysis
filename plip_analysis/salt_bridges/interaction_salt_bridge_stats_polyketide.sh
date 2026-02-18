#!/bin/bash

# --------------------------------------
# User configuration
# --------------------------------------
root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/polyketide"
output_dir="polyketide_saltbridge_stats"
mkdir -p "$output_dir"

# Output files
out_restype="$output_dir/RESTYPE_LIG.txt"
out_dist="$output_dir/DIST.txt"
out_group="$output_dir/LIG_GROUP.txt"

# Clear old content
> "$out_restype"
> "$out_dist"
> "$out_group"

# --------------------------------------
# Process each report.txt
# --------------------------------------
find "$root_dir" -type f -name "report.txt" | while read -r report_path; do

    awk -v f_restype="$out_restype" \
        -v f_dist="$out_dist" \
        -v f_group="$out_group" '
        BEGIN { in_sb = 0 }

        # Enter Salt Bridges section
        /^\*\*Salt Bridges\*\*/ { in_sb = 1; next }

        # Leave Salt Bridges section if another **section** starts
        /^\*\*/ && !/^\*\*Salt Bridges\*\*/ { in_sb = 0 }

        # Skip header row
        in_sb && /^\| RESNR/ { next }

        # Parse each data row
        in_sb && /^\|/ {
            # Remove pipes and compress whitespace
            gsub(/\|/, "", $0)
            gsub(/[ \t]+/, " ", $0)

            split($0, f, " ")

            # Column mapping from PLIP table:
            #  6: RESTYPE_LIG
            #  8: DIST
            # 10: LIG_GROUP

            restype_lig = f[6]
            dist        = f[8]
            lig_group   = f[10]

            if (restype_lig != "") print restype_lig >> f_restype
            if (dist ~ /^[0-9.]+$/) print dist >> f_dist
            if (lig_group != "") print lig_group >> f_group
        }
    ' "$report_path"

done

echo "Salt bridge features extracted to:"
echo "  RESTYPE_LIG -> $out_restype"
echo "  DIST        -> $out_dist"
echo "  LIG_GROUP   -> $out_group"
