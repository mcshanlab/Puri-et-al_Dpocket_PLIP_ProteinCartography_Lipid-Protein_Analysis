#!/bin/bash

root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/Fatty_acyl"
temp_dir="fattyacyl_pistacking_stats"
mkdir -p "$temp_dir"

# Output files
out_centdist="$temp_dir/CENTDIST.txt"
out_angle="$temp_dir/ANGLE.txt"
out_offset="$temp_dir/OFFSET.txt"
out_type="$temp_dir/TYPE.txt"

# Clear old content
> "$out_centdist"
> "$out_angle"
> "$out_offset"
> "$out_type"

find "$root_dir" -type f -name "report.txt" | while read -r report_path; do

    awk -v f_cd="$out_centdist" \
        -v f_ang="$out_angle" \
        -v f_off="$out_offset" \
        -v f_typ="$out_type" '
        
        # trim leading/trailing whitespace
        function trim(s) {
            gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s)
            return s
        }

        BEGIN { in_pi = 0 }

        # Enter pi-stacking section (case-insensitive)
        /^\*\*[Pp]i-Stacking\*\*/ { in_pi = 1; next }

        # Leave when another section begins
        /^\*\*/ && !/^\*\*[Pp]i-Stacking\*\*/ { in_pi = 0 }

        # Process only table rows
        in_pi && /^\|/ {

            # Split on "|" (keep exact columns)
            n = split($0, parts, "|")

            # Need at least 12 pieces to contain required columns
            if (n >= 12) {

                resnr = trim(parts[2])

                # Skip header and separators; only numeric RESNR is data
                if (resnr ~ /^[0-9]+$/) {

                    centdist = trim(parts[9])
                    angle    = trim(parts[10])
                    offset   = trim(parts[11])
                    typecol  = trim(parts[12])

                    if (centdist ~ /^[0-9.]+$/) print centdist >> f_cd
                    if (angle    ~ /^[0-9.]+$/) print angle    >> f_ang
                    if (offset   ~ /^[0-9.]+$/) print offset   >> f_off
                    if (typecol  != "")          print typecol >> f_typ
                }
            }
        }
    ' "$report_path"

done

echo "Pi-stacking features extracted to: $temp_dir"
