#!/bin/bash

# ----------------------------------------
# Config
# ----------------------------------------
root_dir="/Volumes/GigiMurin/plip/sphingo_lipids"
output_file="Z1T_lipid_atoms_total_clean.txt"
temp_file="temp_all_atoms.tmp"
temp_pdb_file="temp_atom_names.tmp"

# Clean previous temp files
> "$temp_file"
> "$temp_pdb_file"
> "$output_file"

# ----------------------------------------
# Traverse all report.txt files in Z1T folders
# ----------------------------------------
find "$root_dir" -type f -name "report.txt" | while read -r report_path; do
    folder_name=$(basename "$(dirname "$report_path")")
    if [[ "$folder_name" =~ Z1T[0-9]+$ ]]; then

        # ---------------------------
        # Parse PDB to get Z1T atom index -> atom name + element
        # ---------------------------
        pdb_file=$(find "$(dirname "$report_path")" -maxdepth 1 -name "*.pdb" | head -n 1)

        if [[ -f "$pdb_file" ]]; then
            awk '
            $1=="HETATM" && $4=="Z1T" {
                idx = substr($0,7,5)+0
                atomname = substr($0,13,4)
                gsub(/ /,"",atomname)
                element = substr($0,77,2)
                gsub(/ /,"",element)
                if(atomname!="" && element!="")
                    print idx, atomname, element
            }' "$pdb_file" >> "$temp_pdb_file"
        fi

        # ---------------------------
        # Parse PLIP report for atom indices
        # ---------------------------
        awk '
            /^\*\*/ { in_table=0; next }

            /^\|/ && (index($0,"LIGCARBONIDX") || index($0,"LIG_IDX_LIST") || index($0,"DONORIDX") || index($0,"ACCEPTORIDX")) {
                in_table=1
                header_line=$0
                gsub(/\|/,"",header_line)
                split(header_line, headers, " ")
                delete col
                for(i=1;i<=length(headers);i++) col[headers[i]]=i
                next
            }

            in_table && /^\|/ {
                gsub(/\|/,"",$0)
                split($0, fields, " ")

                if("LIGCARBONIDX" in col && fields[col["LIGCARBONIDX"]] ~ /^[0-9,]+$/){
                    split(fields[col["LIGCARBONIDX"]], atoms, ",")
                    for(a in atoms) print atoms[a] >> "'"$temp_file"'"
                }

                if("DONORIDX" in col && fields[col["DONORIDX"]] ~ /^[0-9]+$/){
                    print fields[col["DONORIDX"]] >> "'"$temp_file"'"
                }

                if("ACCEPTORIDX" in col && fields[col["ACCEPTORIDX"]] ~ /^[0-9]+$/){
                    print fields[col["ACCEPTORIDX"]] >> "'"$temp_file"'"
                }

                if("LIG_IDX_LIST" in col && fields[col["LIG_IDX_LIST"]] ~ /^[0-9,]+$/){
                    split(fields[col["LIG_IDX_LIST"]], atoms, ",")
                    for(a in atoms) print atoms[a] >> "'"$temp_file"'"
                }
            }
        ' "$report_path"

    fi
done

# ----------------------------------------
# Map atom indices -> name + element
# Initialize ALL atoms with count = 0
# Then increment when seen
# ----------------------------------------
awk '
NR==FNR {
    # Build index -> atomname,element map
    map[$1]=$2","$3

    # Initialize count to 0 for every atom
    key=$2","$3
    if(!(key in count)){
        count[key]=0
    }
    next
}

{
    idx=$1
    if(idx in map){
        split(map[idx],a,",")
        key=a[1]","a[2]
        count[key]++
    }
}

END{
    for(k in count){
        split(k,b,",")
        printf "%s %d %s\n", b[1], count[k], b[2]
    }
}
' "$temp_pdb_file" "$temp_file" | sort -k2,2nr > "$output_file"

# ----------------------------------------
# Clean up
# ----------------------------------------
rm "$temp_file" "$temp_pdb_file"

echo "Consolidated Z1T atom counts saved in: $output_file"
