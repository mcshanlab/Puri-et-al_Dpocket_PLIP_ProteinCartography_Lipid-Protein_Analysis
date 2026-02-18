#!/bin/bash

# ----------------------------------------
# Config
# ----------------------------------------
root_dir="/Volumes/GigiMurin/plip/sphingo_lipids"
output_file="16C_lipid_atoms_total_clean.txt"
temp_file="temp_all_atoms.tmp"
temp_pdb_file="temp_atom_names.tmp"

# Clean previous temp files
> "$temp_file"
> "$temp_pdb_file"
> "$output_file"

# ----------------------------------------
# Traverse all report.txt files in 16C folders
# ----------------------------------------
find "$root_dir" -type f -name "report.txt" | while read -r report_path; do
    folder_name=$(basename "$(dirname "$report_path")")
    if [[ "$folder_name" =~ 16C[0-9]+$ ]]; then

        # ---------------------------
        # Parse PDB to get 16C atom index -> atom name + element
        # ---------------------------
        pdb_file=$(find "$(dirname "$report_path")" -maxdepth 1 -name "*.pdb" | head -n 1)
        if [[ -f "$pdb_file" ]]; then
            awk '{
                # Only HETATM lines for 16C ligand
                if($1=="HETATM" && $4=="16C"){
                    idx = substr($0,7,5)+0
                    atomname = substr($0,13,4)
                    gsub(/ /,"",atomname)
                    element = substr($0,77,2)
                    gsub(/ /,"",element)
                    if(atomname!="" && element!="") print idx, atomname, element
                }
            }' "$pdb_file" >> "$temp_pdb_file"
        fi

        # ---------------------------
        # Parse PLIP report for atom indices (all interaction types)
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

            /^\*\*/ { in_table=0 }

            in_table && /^\|/ {
                gsub(/\|/,"",$0)
                split($0, fields, " ")

                # LIGCARBONIDX
                if("LIGCARBONIDX" in col && fields[col["LIGCARBONIDX"]] ~ /^[0-9,]+$/){
                    split(fields[col["LIGCARBONIDX"]], atoms, ",")
                    for(a in atoms) print atoms[a] >> "'"$temp_file"'"
                }

                # DONORIDX
                if("DONORIDX" in col && fields[col["DONORIDX"]] ~ /^[0-9]+$/){
                    print fields[col["DONORIDX"]] >> "'"$temp_file"'"
                }

                # ACCEPTORIDX
                if("ACCEPTORIDX" in col && fields[col["ACCEPTORIDX"]] ~ /^[0-9]+$/){
                    print fields[col["ACCEPTORIDX"]] >> "'"$temp_file"'"
                }

                # LIG_IDX_LIST
                if("LIG_IDX_LIST" in col && fields[col["LIG_IDX_LIST"]] ~ /^[0-9,]+$/){
                    split(fields[col["LIG_IDX_LIST"]], atoms, ",")
                    for(a in atoms) print atoms[a] >> "'"$temp_file"'"
                }
            }
        ' "$report_path"

    fi
done

# ----------------------------------------
# Map atom indices -> name + element, sum counts by atom name
# ----------------------------------------
awk 'NR==FNR { map[$1]=$2","$3; next }
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
     }' "$temp_pdb_file" "$temp_file" | sort -k2,2nr > "$output_file"

# ----------------------------------------
# Clean up
# ----------------------------------------
rm "$temp_file" "$temp_pdb_file"

echo "Consolidated 16C atom counts saved in: $output_file"
