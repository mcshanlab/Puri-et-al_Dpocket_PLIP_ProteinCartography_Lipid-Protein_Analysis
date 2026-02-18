#!/bin/bash

root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/sphingo_lipids"
temp_dir="sphingolipid_halogenbond_stats"
mkdir -p "$temp_dir"

# Output files
out_sidechain="$temp_dir/sidechain.txt"
out_dist="$temp_dir/dist.txt"
out_don_ang="$temp_dir/donor_angle.txt"
out_acc_ang="$temp_dir/acceptor_angle.txt"
out_donor_type="$temp_dir/donor_type.txt"
out_acceptor_type="$temp_dir/acceptor_type.txt"

# Clear old contents
> "$out_sidechain"
> "$out_dist"
> "$out_don_ang"
> "$out_acc_ang"
> "$out_donor_type"
> "$out_acceptor_type"

# Process all PLIP reports
find "$root_dir" -type f -name "report.txt" | while read -r report_path; do

    awk \
        -v f_side="$out_sidechain" \
        -v f_dist="$out_dist" \
        -v f_dang="$out_don_ang" \
        -v f_aang="$out_acc_ang" \
        -v f_dt="$out_donor_type" \
        -v f_at="$out_acceptor_type" '

        BEGIN { in_halo = 0 }

        # Enter halogen bond section
        /^\*\*Halogen Bonds\*\*/ { in_halo = 1; next }

        # Leave halogen bond section when a new **section** begins
        /^\*\*/ && !/^\*\*Halogen Bonds\*\*/ { in_halo = 0 }

        # Skip header
        in_halo && /^\| RESNR/ { next }

        # Extract table rows
        in_halo && /^\|/ {
            line = $0

            # remove pipes, collapse whitespace
            gsub(/\|/, "", line)
            gsub(/[ \t]+/, " ", line)

            split(line, f, " ")

            # Column map for Halogen Bond table:
            #   7: SIDECHAIN
            #   8: DIST
            #   9: DON_ANGLE
            #  10: ACC_ANGLE
            #  12: DONORTYPE
            #  14: ACCEPTORTYPE

            side = f[7]
            dist = f[8]
            dang = f[9]
            aang = f[10]
            dtyp = f[12]
            atyp = f[14]

            if (side != "") print side >> f_side
            if (dist ~ /^[0-9.]+$/) print dist >> f_dist
            if (dang ~ /^[0-9.]+$/) print dang >> f_dang
            if (aang ~ /^[0-9.]+$/) print aang >> f_aang
            if (dtyp != "") print dtyp >> f_dt
            if (atyp != "") print atyp >> f_at
        }

    ' "$report_path"

done

echo "Halogen bond features extracted to: $temp_dir"
