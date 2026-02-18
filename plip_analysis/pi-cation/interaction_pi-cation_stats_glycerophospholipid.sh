#!/bin/bash

root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/Glycerophospholipid"
output_dir="glycerophospholipid_pication_stats"
mkdir -p "$output_dir"

# Output files
out_dist="$output_dir/DIST.txt"
out_offset="$output_dir/OFFSET.txt"
out_protcharged="$output_dir/PROTCHARGED.txt"
out_liggroup="$output_dir/LIG_GROUP.txt"

# Clear old content
> "$out_dist"
> "$out_offset"
> "$out_protcharged"
> "$out_liggroup"

find "$root_dir" -type f -name "report.txt" | while read -r report_path; do

    awk -v f_dist="$out_dist" \
        -v f_off="$out_offset" \
        -v f_pc="$out_protcharged" \
        -v f_lg="$out_liggroup" '
        
        BEGIN { in_pc = 0 }

        # Enter pi-Cation section
        /^\*\*pi-Cation/ { in_pc = 1; next }

        # Leave section if another **section** starts
        /^\*\*/ && !/^\*\*pi-Cation/ { in_pc = 0 }

        # Skip header row
        in_pc && /^\| RESNR/ { next }

        # Process data rows
        in_pc && /^\|/ {
            # Remove pipes and compress spaces
            gsub(/\|/, "", $0)
            gsub(/[ \t]+/, " ", $0)
            $0 = trim($0)

            split($0, f, " ")

            # Column mapping from PLIP pi-Cation table:
            #  8: DIST
            #  9: OFFSET
            # 10: PROTCHARGED
            # 11: LIG_GROUP
            dist        = f[8]
            offset      = f[9]
            protcharged = f[10]
            liggroup    = f[11]

            if (dist ~ /^[0-9.]+$/)        print dist        >> f_dist
            if (offset ~ /^[0-9.]+$/)      print offset      >> f_off
            if (protcharged != "")          print protcharged >> f_pc
            if (liggroup != "")             print liggroup    >> f_lg
        }

        # Trim function
        function trim(s) { gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s); return s }

    ' "$report_path"

done

echo "π-Cation interaction features extracted to: $output_dir"
echo "  DIST        -> $out_dist"
echo "  OFFSET      -> $out_offset"
echo "  PROTCHARGED -> $out_protcharged"
echo "  LIG_GROUP   -> $out_liggroup"
