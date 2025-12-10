#!/bin/bash

# deal wiht all species in batch

cd /rhome/yyen008/bigdata/gen220/Project_Una_2025

# define species array
declare -a species=(
    "Phocaeicola_dorei_HM719"
    "Phocaeicola_vulgatus_ATCC"
    "Bacteroides_thetaiotaomicron_VPI-5482"
    "Bacteroides_thetaiotaomicron_jmh43"
    "Bacteroides_intestinalis_APC919"
    "Bacteroides_uniformis_ATCC"
    "Parabacteroides_distasonis_DA491"
    "Alistipes_shahii_DA73"
    "Escherichia_coli_MG1655"
    "Escherichia_coli_BL21"
)

echo "=========================================="
echo "start dealing with all species"
echo "=========================================="
echo ""

for species in "${species[@]}"; do
    echo "deal with: $species"
    python3 filter_ortholog_sequences.py "$species" "genome_seq" "core_gene"
    echo ""
done

echo "=========================================="
echo "completed dealing with all species"
echo "=========================================="

# summarize output
echo ""
echo "output summary:"
echo "---"
cd core_gene
for file in *_core_protein.faa; do
    species_name=$(echo $file | sed 's/_core_protein.faa//')
    core_prot=$(grep -c "^>" "${species_name}_core_protein.faa" 2>/dev/null || echo "0")
    other_prot=$(grep -c "^>" "${species_name}_other_protein.faa" 2>/dev/null || echo "0")
    core_cds=$(grep -c "^>" "${species_name}_core_CDS.fna" 2>/dev/null || echo "0")
    other_cds=$(grep -c "^>" "${species_name}_other_CDS.fna" 2>/dev/null || echo "0")
    
    printf "%-40s | Core: %5d | Other: %5d | Protein: %5d | %5d\n" \
        "$species_name" "$core_cds" "$other_cds" "$core_prot" "$other_prot"
done

echo ""
echo "attention："
echo "- CDS column = Coding DNA Sequence (including pseudo genes)"
echo "- Protein column = Functional proteins (not including pseudo genes)"
echo "- if the number of CDS and Protein are differnet, should caused by pseudo genes"
