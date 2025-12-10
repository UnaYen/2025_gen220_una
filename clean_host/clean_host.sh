#!/usr/bin/env bash
set -euo pipefail
module load ncbi-blast/2.11.0+
module load samtools/1.14
module load bedtools/2.30.0
module load seqkit/2.4.0
#############################################
# 找出 phage 在 host genome 中的整段整合區域，
# 把那一段從 genome 裡去掉，然後輸出：
#   1) Host_clean_genomic.fna  (genome 去掉 phage 區段，含單一 header 和 merged sequence)
#   2) Host_clean_CDS.fna      (只含 host 的 CDS nucleotide)
#   3) Host_clean.fnn          (只含 host 的 CDS nucleotide - 原始格式)
#   4) Host_clean.faa          (只含 host 的 protein)
#
# 使用方法：
#   bash clean_host.sh -input HOST_PREFIX PHAGE_FNA -output all
#
# 其中 HOST_PREFIX 會自動拼成：
#   HOST_PREFIX_genomic.fna
#   HOST_PREFIX_genomic.gff
#   HOST_PREFIX_protein.faa
#############################################

usage() {
    echo "Usage:"
    echo "  bash $0 -input HOST_FNA PHAGE_FNA [-minlen LENGTH] -output all"
    echo
    echo "Example:"
    echo "  bash $0 -input ../genome_seq/Phocaeicola_vulgatus_ATCC_genomic.fna ../genome_seq/Pv_PV01_genomic.fna -minlen 3000 -output all"
    echo
    echo "Parameters:"
    echo "  -input       HOST_FNA and PHAGE_FNA files (required)"
    echo "  -minlen      Minimum alignment length in bp (default: 3000)"
    echo "  -output      Output mode (default: all)"
    echo
    exit 1
}

if [[ $# -lt 1 ]]; then
    usage
fi

HOST_PREFIX=""
HOST_FNA=""
HOST_GFF=""
HOST_FAA=""
HOST_CDS_FNA=""
PHAGE_FNA=""
OUTPUT_MODE="all"
MIN_ALIGN_LENGTH=3000

# 簡單的參數 parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        -input)
            # -input 後面應該接兩個東西：HOST_FNA PHAGE_FNA
            if [[ $# -lt 3 ]]; then
                echo "ERROR: -input 需要兩個參數：HOST_FNA PHAGE_FNA" >&2
                usage
            fi
            HOST_FNA="$2"
            PHAGE_FNA="$3"
            shift 3
            ;;
        -minlen)
            # -minlen 後面應該接一個數字
            if [[ $# -lt 2 ]]; then
                echo "ERROR: -minlen 需要一個參數（最小比對長度）" >&2
                usage
            fi
            MIN_ALIGN_LENGTH="$2"
            shift 2
            ;;
        -output)
            if [[ $# -lt 2 ]]; then
                echo "ERROR: -output 需要一個參數（目前只接受 all）" >&2
                usage
            fi
            OUTPUT_MODE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option or extra argument: $1" >&2
            usage
            ;;
    esac
done

if [[ -z "${HOST_FNA}" || -z "${PHAGE_FNA}" ]]; then
    echo "ERROR: 必須指定 -input HOST_FNA PHAGE_FNA" >&2
    usage
fi

if [[ "${OUTPUT_MODE}" != "all" ]]; then
    echo "ERROR: 目前只支援 -output all" >&2
    exit 1
fi

#############################################
# 從 HOST_FNA 推導 HOST_PREFIX 和其他檔案名稱
#############################################

# 獲得 directory 和 basename
HOST_DIR=$(dirname "${HOST_FNA}")
HOST_BASENAME=$(basename "${HOST_FNA}")

# 去掉 .fna 或 _genomic.fna 後綴，得到 prefix
if [[ "${HOST_BASENAME}" == *"_genomic.fna" ]]; then
    HOST_PREFIX="${HOST_BASENAME%_genomic.fna}"
else
    HOST_PREFIX="${HOST_BASENAME%.fna}"
fi

# 在同一個 directory 裡找對應的 GFF 和 FAA 檔案
HOST_GFF="${HOST_DIR}/${HOST_PREFIX}_genomic.gff"
HOST_FAA="${HOST_DIR}/${HOST_PREFIX}_protein.faa"

# 自動搜尋 CDS 檔案 - 尋找 *cds_from_genomic.fna 檔案
# 首先嘗試用 HOST_PREFIX 找對應的 CDS 檔案
if [[ -f "${HOST_DIR}/${HOST_PREFIX}_cds_from_genomic.fna" ]]; then
    HOST_CDS_FNA="${HOST_DIR}/${HOST_PREFIX}_cds_from_genomic.fna"
else
    # 如果找不到，嘗試在同一個 directory 裡找任何 *cds_from_genomic.fna 檔案
    CDS_FILES=$(find "${HOST_DIR}" -maxdepth 1 -name "*cds_from_genomic.fna" 2>/dev/null | head -1)
    if [[ -n "${CDS_FILES}" ]]; then
        HOST_CDS_FNA="${CDS_FILES}"
    fi
fi

#############################################
# 檢查 host 的各種檔案是否存在
#############################################

if [[ ! -f "${HOST_FNA}" ]]; then
    echo "ERROR: 找不到 host genome 檔案: ${HOST_FNA}" >&2
    exit 1
fi
if [[ ! -f "${HOST_GFF}" ]]; then
    echo "ERROR: 找不到 host GFF 註解檔: ${HOST_GFF}" >&2
    exit 1
fi
if [[ ! -f "${HOST_FAA}" ]]; then
    echo "ERROR: 找不到 host protein 檔案: ${HOST_FAA}" >&2
    exit 1
fi
if [[ ! -f "${PHAGE_FNA}" ]]; then
    echo "ERROR: 找不到 phage genome 檔案: ${PHAGE_FNA}" >&2
    exit 1
fi

# 檢查 CDS 檔案是否存在
if [[ -n "${HOST_CDS_FNA}" ]]; then
    if [[ ! -f "${HOST_CDS_FNA}" ]]; then
        echo "ERROR: 找不到 host CDS 檔案: ${HOST_CDS_FNA}" >&2
        exit 1
    fi
    echo "==> 找到 CDS 檔案：${HOST_CDS_FNA}"
else
    echo "==> 未找到 CDS 檔案，將使用 bedtools 從 host genome 抽取"
fi
if [[ ! -f "${PHAGE_FNA}" ]]; then
    echo "ERROR: 找不到 phage genome 檔案: ${PHAGE_FNA}" >&2
    exit 1
fi

#############################################
# 檢查必要程式是否存在
#############################################

for cmd in makeblastdb blastn samtools bedtools seqkit; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: 找不到指令：$cmd ；請先在環境中安裝或 module load" >&2
        exit 1
    fi
done

#############################################
# 設定輸出檔案名稱
#############################################

OUT_PREFIX="${HOST_PREFIX}_Host_clean"

OUT_FNA_GENOMIC="${OUT_PREFIX}_genomic.fna"
OUT_FNA_CDS="${OUT_PREFIX}_CDS.fna"
OUT_FNN="${OUT_PREFIX}.fnn"
OUT_FAA="${OUT_PREFIX}.faa"

BLAST_DB="${HOST_PREFIX}_BLASTDB"
BLAST_OUT="${HOST_PREFIX}_phage_vs_host.tab"
PHAGE_BED="phage_region.bed"
HOST_GENOME_SIZES="Host.genome"
HOST_CDS_BED="Host_CDS.bed"
HOST_CDS_CLEAN_BED="Host_clean_CDS.bed"
HOST_CLEAN_REGIONS_BED="Host_clean_regions.bed"
HOST_CLEAN_IDS="Host_clean_ids.txt"

#############################################
# 第一步：建立 host 的 BLAST 資料庫
#############################################

echo "==> 建立 BLAST 資料庫：${HOST_FNA}"
makeblastdb -in "${HOST_FNA}" -dbtype nucl -out "${BLAST_DB}" >/dev/null

#############################################
# 第二步：用 phage genome 去打 BLAST
#############################################

echo "==> BLAST phage (${PHAGE_FNA}) 對 host (${HOST_FNA})"
blastn \
    -query "${PHAGE_FNA}" \
    -db "${BLAST_DB}" \
    -out "${BLAST_OUT}" \
    -outfmt 6 \
    -max_target_seqs 50 \
    -max_hsps 50

if [[ ! -s "${BLAST_OUT}" ]]; then
    echo "ERROR: BLAST 沒有任何 hit，無法自動界定 phage 區域。" >&2
    exit 1
fi

#############################################
# 第三步：從 BLAST 結果推算 phage 在 host 上的大區段
#
# 策略：
#   1) 只考慮比對長度 > 3000 bp 的 hit
#   2) 對每個 host contig，計算符合條件的 phage 對上去的總長度 total_hit_len
#   3) 找出 total_hit_len 最大的那條 contig
#   4) 在那條 contig 上，取所有符合條件的 hit 的最小座標做 start，
#      最大座標做 end，當成 prophage 的整體範圍
#############################################

echo "==> 從 BLAST 結果自動推定 phage 整合區域 (最小長度: ${MIN_ALIGN_LENGTH} bp)"

awk -v min_len="${MIN_ALIGN_LENGTH}" '
{
    len = $4        # alignment length
    
    # 只考慮長度 > min_len 的 hit
    if (len <= min_len) {
        next
    }
    
    cont = $2       # sseqid = host contig
    # sstart / send 可能大小顛倒，先整理成 s <= e
    s = ($9 < $10 ? $9 : $10)
    e = ($9 > $10 ? $9 : $10)

    total_len[cont] += len

    if (!(cont in min_s) || s < min_s[cont]) {
        min_s[cont] = s
    }
    if (!(cont in max_e) || e > max_e[cont]) {
        max_e[cont] = e
    }
}
END {
    best_contig = ""
    best_len = 0
    for (c in total_len) {
        if (total_len[c] > best_len) {
            best_len = total_len[c]
            best_contig = c
        }
    }
    if (best_contig == "") {
        exit 1
    }
    # BED 是 0-based，因此 start 要減 1，但不能小於 0
    start = min_s[best_contig] - 1
    if (start < 0) start = 0
    end = max_e[best_contig]

    # 輸出：contig, start, end
    printf("%s\t%d\t%d\n", best_contig, start, end)
}
' "${BLAST_OUT}" > "${PHAGE_BED}"

if [[ ! -s "${PHAGE_BED}" ]]; then
    echo "ERROR: 無法從 BLAST 結果推定 phage 區域 (沒有比對長度 > ${MIN_ALIGN_LENGTH} bp 的 hit)。" >&2
    exit 1
fi

echo "    推定 phage 區域 (BED)："
cat "${PHAGE_BED}"

#############################################
# 第四步：建立 genome 長度檔案，算出非-phage 區域
#############################################

echo "==> 建立 genome index 並計算非-phage 區域"

samtools faidx "${HOST_FNA}"
cut -f1,2 "${HOST_FNA}.fai" > "${HOST_GENOME_SIZES}"

bedtools complement -i "${PHAGE_BED}" -g "${HOST_GENOME_SIZES}" > "${HOST_CLEAN_REGIONS_BED}"

#############################################
# 第五步：從原始 genome 抓出非-phage 的序列 -> Host_clean_genomic.fna
# 使用 bedtools getfasta，然後合併所有序列成單一 header 和單一序列
#############################################

echo "==> 產生臨時 ${OUT_FNA_GENOMIC}（去掉 phage 區域的 host genome）"

bedtools getfasta \
    -fi "${HOST_FNA}" \
    -bed "${HOST_CLEAN_REGIONS_BED}" \
    -fo "${OUT_FNA_GENOMIC}.tmp"

# 合併所有序列為單一 header 和單一行序列
echo "==> 合併序列成單一 header 和單一序列"
awk '
BEGIN {
    seq = ""
    header = ""
}
/^>/ {
    if (header == "") {
        # 第一個 header，取第一個序列的 ID 去掉冒號後面的部分
        header = substr($0, 2)
        # 只保留第一個 ID（冒號前面）
        idx = index(header, ":")
        if (idx > 0) {
            header = substr(header, 1, idx - 1)
        }
    }
    # 跳過後續的 header，只取第一個
    next
}
{
    # 累積序列
    seq = seq $0
}
END {
    if (header != "" && seq != "") {
        printf(">%s\n%s\n", header, seq)
    }
}
' "${OUT_FNA_GENOMIC}.tmp" > "${OUT_FNA_GENOMIC}"

rm -f "${OUT_FNA_GENOMIC}.tmp"

#############################################
# 第六步：從 GFF 拿出所有 CDS -> BED
#############################################

echo "==> 從 GFF 轉成 CDS BED"

# 檢查是否是 Bacteroides 格式（source 欄位為 "Protein Homology"）
BACTEROIDES_FORMAT=0
head -1000 "${HOST_GFF}" > /tmp/gff_check.txt || true
if grep "Protein.*Homology.*CDS" /tmp/gff_check.txt > /dev/null 2>&1; then
    BACTEROIDES_FORMAT=1
fi
rm -f /tmp/gff_check.txt

if [[ ${BACTEROIDES_FORMAT} -eq 1 ]]; then
    echo "    偵測到 Bacteroides 格式的 GFF（source 為 'Protein Homology'）"
    
    # 使用特殊的轉換邏輯
    # Bacteroides GFF 格式：Col1=seqid, Col2=Protein, Col3=Homology, Col4=CDS, Col5=start, Col6=end, Col7=., Col8=strand, Col9=phase, Col10=attributes
    awk '
    $2 == "Protein" && $3 == "Homology" && $4 == "CDS" {
        start = $5 - 1  # GFF 是 1-based，BED 需要 0-based
        if (start < 0) start = 0
        end = $6
        strand = $8
        attr = $10
        
        printf("%s\t%d\t%d\t%s\t0\t%s\n", $1, start, end, attr, strand)
    }
    ' "${HOST_GFF}" > "${HOST_CDS_BED}"
else
    # 使用標準 GFF 轉換
    awk '
    $3 == "CDS" {
        # GFF: col4, col5 是 1-based，BED 需要 0-based，所以下限要 -1
        start = $4 - 1
        if (start < 0) start = 0
        printf("%s\t%d\t%s\t%s\t0\t%s\n", $1, start, $5, $9, $7)
    }
    ' "${HOST_GFF}" > "${HOST_CDS_BED}"
fi

#############################################
# 第七步：移除 overlap 到 phage 區域的 CDS
#############################################

echo "==> 移除落在 phage 區域裡的 CDS"

bedtools subtract \
    -a "${HOST_CDS_BED}" \
    -b "${PHAGE_BED}" \
    > "${HOST_CDS_CLEAN_BED}"

#############################################
# 第八步：產生乾淨的 CDS 序列 -> Host_clean_CDS.fna
#############################################

if [[ -n "${HOST_CDS_FNA}" ]]; then
    # 如果提供了 CDS 檔案，則從 CDS 檔案中篩選出不在 phage 區域的 CDS
    echo "==> 從 CDS 檔案篩選出不在 phage 區域的 CDS -> ${OUT_FNA_CDS}"
    
    # 從 HOST_CDS_CLEAN_BED 中提取 CDS ID（去掉 cds- 前綴），用於篩選
    grep -o "ID=[^;]*" "${HOST_CDS_CLEAN_BED}" | sed "s/ID=//" | sed "s/^cds-//" | sort -u > "${HOST_CLEAN_IDS}"
    
    # 用 awk 從 CDS 檔案中篩選出在 ID 清單中的序列
    # CDS 檔案的 header 包含 [protein_id=XXX] 格式
    awk '
    BEGIN {
        # 讀取 ID 清單
        while ((getline id < "'"${HOST_CLEAN_IDS}"'") > 0) {
            ids[id] = 1
        }
        close("'"${HOST_CLEAN_IDS}"'")
    }
    /^>/ {
        # 從 header 中提取 protein_id=XXX
        if (match($0, /protein_id=([^ \]]+)/, arr)) {
            protein_id = arr[1]
            if (protein_id in ids) {
                print_seq = 1
            } else {
                print_seq = 0
            }
        }
        if (print_seq) print $0
        next
    }
    {
        if (print_seq) print $0
    }
    ' "${HOST_CDS_FNA}" > "${OUT_FNA_CDS}"
else
    # 如果沒有提供 CDS 檔案，則從 host genome 中用 bedtools 抓出
    echo "==> 從 host genome 抓出乾淨的 CDS 序列 -> ${OUT_FNA_CDS}"
    
    bedtools getfasta \
        -fi "${HOST_FNA}" \
        -bed "${HOST_CDS_CLEAN_BED}" \
        -s \
        -name \
        > "${OUT_FNA_CDS}"
fi

#############################################
# 第八步B：產生原始格式的 CDS 序列 -> Host_clean.fnn
#############################################

echo "==> 產生 ${OUT_FNN}（只含 host CDS 的 nucleotide - 原始格式）"

bedtools getfasta \
    -fi "${HOST_FNA}" \
    -bed "${HOST_CDS_CLEAN_BED}" \
    -s \
    -name \
    > "${OUT_FNN}"

#############################################
# 第九步：根據 CDS ID 清出乾淨的 protein -> Host_clean.faa
#############################################

echo "==> 產生 ${OUT_FAA}（只含 host 的 protein）"

# 從 BED 的第 4 欄把 ID=xxx 抓出來，再去掉 ID=，並去掉 cds- 前綴
# 如果之前已經生成過（CDS 檔案情況），跳過此步驟
if [[ ! -f "${HOST_CLEAN_IDS}" ]]; then
    grep -o "ID=[^;]*" "${HOST_CDS_CLEAN_BED}" | sed "s/ID=//" | sed "s/^cds-//" | sort -u > "${HOST_CLEAN_IDS}"
fi

# 根據這份 ID 清出 protein
seqkit grep -f "${HOST_CLEAN_IDS}" "${HOST_FAA}" > "${OUT_FAA}"

echo "==> 完成！輸出檔案："
echo "    ${OUT_FNA_GENOMIC}"
echo "    ${OUT_FNA_CDS}"
echo "    ${OUT_FNN}"
echo "    ${OUT_FAA}"
