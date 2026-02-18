#!/bin/bash

root_dir="/Volumes/8TB_McShan_Drive/nikki/project_1/plip/prenol_lipid"
output_dir="prenol_metalcomplex_stats"
mkdir -p "$output_dir"

# Output files
out_metaltype="$output_dir/METAL_TYPE.txt"
out_targettype="$output_dir/TARGET_TYPE.txt"
out_coord="$output_dir/COORDINATION.txt"
out_dist="$output_dir/DIST.txt"
out_location="$output_dir/LOCATION.txt"
out_geometry="$output_dir/GEOMETRY.txt"
out_restype_lig="$output_dir/RESTYPE_LIG.txt"

# Clear old content
> "$out_metaltype"
> "$out_targettype"
> "$out_coord"
> "$out_dist"
> "$out_location"
> "$out_geometry"
> "$out_restype_lig"

find "$root_dir" -type f -name "report.txt" | while read -r report_path; do

    awk \
        -v f_mt="$out_metaltype" \
        -v f_tt="$out_targettype" \
        -v f_co="$out_coord" \
        -v f_di="$out_dist" \
        -v f_lo="$out_location" \
        -v f_ge="$out_geometry" \
        -v f_rl="$out_restype_lig" '

        BEGIN { in_mc = 0 }

        # Enter Metal Complexes section
        /^\*\*Metal Complexes\*\*/ { in_mc = 1; next }

        # Leave when a new section starts
        /^\*\*/ && !/^\*\*Metal Complexes\*\*/ { in_mc = 0 }

        # Skip table header
        in_mc && /^\| RESNR/ { next }

        # Parse metal complex table rows
        in_mc && /^\|/ {
            line = $0

            # Remove pipes and normalize whitespace
            gsub(/\|/, "", line)
            gsub(/[ \t]+/, " ", line)

            split(line, f, " ")

            # PLIP Metal Complex column mapping (after cleanup):
            #  5 : RESTYPE_LIG
            #  8 : METAL_TYPE
            # 10 : TARGET_TYPE
            # 11 : COORDINATION
            # 12 : DIST
            # 13 : LOCATION
            # 15 : GEOMETRY

            restype_lig = f[5]
            metaltype   = f[8]
            targettype  = f[10]
            coord       = f[11]
            dist        = f[12]
            location    = f[13]
            geometry    = f[15]

            if (restype_lig != "") print restype_lig >> f_rl
            if (metaltype != "") print metaltype >> f_mt
            if (targettype != "") print targettype >> f_tt
            if (coord ~ /^[0-9]+$/) print coord >> f_co
            if (dist ~ /^[0-9.]+$/) print dist >> f_di
            if (location != "") print location >> f_lo
            if (geometry != "") print geometry >> f_ge
        }

    ' "$report_path"

done

echo "Metal Complex features extracted to: $output_dir"
echo "  METAL_TYPE    -> $out_metaltype"
echo "  TARGET_TYPE   -> $out_targettype"
echo "  COORDINATION  -> $out_coord"
echo "  DIST          -> $out_dist"
echo "  LOCATION      -> $out_location"
echo "  GEOMETRY      -> $out_geometry"
echo "  RESTYPE_LIG   -> $out_restype_lig"
