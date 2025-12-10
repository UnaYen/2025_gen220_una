#!/usr/bin/env bash
set -euo pipefail

#############################################
# 專門處理 Bacteroides GFF 格式的轉換
# Bacteroides GFF 的 source 欄位是 "Protein Homology"（兩個詞）
# 導致標準的 $3 == "CDS" 失效
#############################################

usage() {
    echo "Usage:"
    echo "  bash $0 <input_gff_file> <output_bed_file>"
    echo
    echo "Example:"
    echo "  bash $0 Bacteroides_thetaiotaomicron_jmh43_genomic.gff Host_CDS.bed"
    echo
    exit 1
}

if [[ $# -lt 2 ]]; then
    usage
fi

INPUT_GFF="$1"
OUTPUT_BED="$2"

if [[ ! -f "${INPUT_GFF}" ]]; then
    echo "ERROR: 找不到輸入檔案 ${INPUT_GFF}" >&2
    exit 1
fi

echo "轉換 Bacteroides GFF 格式..."
echo "輸入: ${INPUT_GFF}"
echo "輸出: ${OUTPUT_BED}"

# Bacteroides GFF 格式：
# Col1: seqid
# Col2-3: "Protein Homology"（兩個詞，所以會被 split 成兩個欄位）
# Col3 (實際): CDS
# Col4 (實際): start
# Col5 (實際): end
# Col6 (實際): score (.)
# Col7 (實際): strand (+/-)
# Col8 (實際): phase (0)
# Col9 (實際): attributes (ID=...)

# 方法：搜尋包含 "Protein Homology" 和 "CDS" 的行，
# 然後從該行中提取 CDS 的位置資訊

awk '
$2 == "Protein" && $3 == "Homology" {
    # 這是 Bacteroides 格式
    # CDS 應該在 $4
    if ($4 == "CDS") {
        start = $5 - 1  # GFF 是 1-based，BED 需要 0-based
        if (start < 0) start = 0
        end = $6
        strand = $8
        attr = $NF
        
        printf("%s\t%d\t%s\t%s\t0\t%s\n", $1, start, end, attr, strand)
    }
}
' "${INPUT_GFF}" > "${OUTPUT_BED}"

line_count=$(wc -l < "${OUTPUT_BED}")
echo "完成！共轉換 ${line_count} 條 CDS 記錄"
