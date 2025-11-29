#!/bin/bash

# 批次處理所有物種的 ortholog 序列分類

cd /rhome/yyen008/bigdata/gen220/Project_Una_2025

# 定義所有物種列表
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
echo "開始批次處理所有物種"
echo "=========================================="
echo ""

for species in "${species[@]}"; do
    echo "處理: $species"
    python3 filter_ortholog_sequences.py "$species" "genome_seq" "core_gene"
    echo ""
done

echo "=========================================="
echo "完成所有物種處理"
echo "=========================================="

# 統計所有輸出檔案
echo ""
echo "輸出檔案統計:"
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
echo "注意："
echo "- CDS 欄位 = Coding DNA Sequence (包含 pseudo genes)"
echo "- Protein 欄位 = Functional proteins (不包含 pseudo genes)"
echo "- 若 CDS 和 Protein 數量不同，差異是因為存在 pseudo genes"
