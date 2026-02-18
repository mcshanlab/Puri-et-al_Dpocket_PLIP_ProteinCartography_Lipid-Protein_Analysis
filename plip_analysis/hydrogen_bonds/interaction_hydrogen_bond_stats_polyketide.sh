#!/bin/bash

root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/polyketide"
temp_dir="polyketide_hbond_stats"
mkdir -p "$temp_dir"

# Output files
out_sidechain="$temp_dir/sidechain.txt"
out_dist_HA="$temp_dir/dist_HA.txt"
out_dist_DA="$temp_dir/dist_DA.txt"
out_angle="$temp_dir/donor_angle.txt"
out_donor_type="$temp_dir/donor_type.txt"
out_acceptor_type="$temp_dir/acceptor_type.txt"

# Clear old content
> "$out_sidechain"
> "$out_dist_HA"
> "$out_dist_DA"
> "$out_angle"
> "$out_donor_type"
> "$out_acceptor_type"

find "$root_dir" -type f -name "report.txt" | while read -r report_path; do
    
    awk -v f_side="$out_sidechain" \
        -v f_ha="$out_dist_HA" \
        -v f_da="$out_dist_DA" \
        -v f_ang="$out_angle" \
        -v f_dt="$out_donor_type" \
        -v f_at="$out_acceptor_type" '
        
        BEGIN { in_hb = 0 }

        # Enter Hydrogen Bond section
        /^\*\*Hydrogen Bonds\*\*/ { in_hb = 1; next }

        # Leave Hydrogen Bond section if another **section** starts
        /^\*\*/ && !/^\*\*Hydrogen Bonds\*\*/ { in_hb = 0 }

        # Skip header row
        in_hb && /^\| RESNR/ { next }

        # Parse each data row
        in_hb && /^\|/ {
            # Remove pipes and compress whitespace
            gsub(/\|/, "", $0)
            gsub(/[ \t]+/, " ", $0)

            split($0, f, " ")

            # Column mapping from PLIP table:
            #  7: SIDECHAIN
            #  8: DIST_H-A
            #  9: DIST_D-A
            # 10: DON_ANGLE
            # 13: DONORTYPE
            # 15: ACCEPTORTYPE

            side = f[7]
            dha  = f[8]
            dda  = f[9]
            ang  = f[10]
            dtyp = f[13]
            atyp = f[15]

            if (side != "") print side >> f_side
            if (dha ~ /^[0-9.]+$/) print dha >> f_ha
            if (dda ~ /^[0-9.]+$/) print dda >> f_da
            if (ang ~ /^[0-9.]+$/) print ang >> f_ang
            if (dtyp != "") print dtyp >> f_dt
            if (atyp != "") print atyp >> f_at
        }
    ' "$report_path"

done

echo "Hydrogen bond features extracted to: $temp_dir"
