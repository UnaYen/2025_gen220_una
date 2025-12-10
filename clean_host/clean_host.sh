#!/usr/bin/env bash
set -euo pipefail
module load ncbi-blast/2.11.0+
module load samtools/1.14
module load bedtools/2.30.0
module load seqkit/2.4.0

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


while [[ $# -gt 0 ]]; do
    case "$1" in
        -input)
    
            if [[ $# -lt 3 ]]; then
                echo "ERROR: -input need two inputfile：HOST_FNA PHAGE_FNA" >&2
                usage
            fi
            HOST_FNA="$2"
            PHAGE_FNA="$3"
            shift 3
            ;;
        -minlen)
            
            if [[ $# -lt 2 ]]; then
                echo "ERROR: -minlen number (minlength of blast results）" >&2
                usage
            fi
            MIN_ALIGN_LENGTH="$2"
            shift 2
            ;;
        -output)
            if [[ $# -lt 2 ]]; then
                echo "ERROR: -output （all）" >&2
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
    echo "ERROR: required -input HOST_FNA PHAGE_FNA" >&2
    usage
fi

if [[ "${OUTPUT_MODE}" != "all" ]]; then
    echo "ERROR: set -output all" >&2
    exit 1
fi

#############################################
# from HOST_FNA to get HOST_PREFIX and other filename
#############################################

#  get directory and basename
HOST_DIR=$(dirname "${HOST_FNA}")
HOST_BASENAME=$(basename "${HOST_FNA}")

# delete .fna or _genomic.fna, get prefix
if [[ "${HOST_BASENAME}" == *"_genomic.fna" ]]; then
    HOST_PREFIX="${HOST_BASENAME%_genomic.fna}"
else
    HOST_PREFIX="${HOST_BASENAME%.fna}"
fi

# find GFF and FAA file in same directory 
HOST_GFF="${HOST_DIR}/${HOST_PREFIX}_genomic.gff"
HOST_FAA="${HOST_DIR}/${HOST_PREFIX}_protein.faa"

# search CDS file - *cds_from_genomic.fna 
if [[ -f "${HOST_DIR}/${HOST_PREFIX}_cds_from_genomic.fna" ]]; then
    HOST_CDS_FNA="${HOST_DIR}/${HOST_PREFIX}_cds_from_genomic.fna"
else
    # if not there, try find any *cds_from_genomic.fna file in same directory
    CDS_FILES=$(find "${HOST_DIR}" -maxdepth 1 -name "*cds_from_genomic.fna" 2>/dev/null | head -1)
    if [[ -n "${CDS_FILES}" ]]; then
        HOST_CDS_FNA="${CDS_FILES}"
    fi
fi

#############################################
# check host file existence 
#############################################

if [[ ! -f "${HOST_FNA}" ]]; then
    echo "ERROR: can't find host genome file: ${HOST_FNA}" >&2
    exit 1
fi
if [[ ! -f "${HOST_GFF}" ]]; then
    echo "ERROR: can't find host GFF file: ${HOST_GFF}" >&2
    exit 1
fi
if [[ ! -f "${HOST_FAA}" ]]; then
    echo "ERROR: can't find host protein dile: ${HOST_FAA}" >&2
    exit 1
fi
if [[ ! -f "${PHAGE_FNA}" ]]; then
    echo "ERROR: can't find phage genome file: ${PHAGE_FNA}" >&2
    exit 1
fi

# check CDS file existence
if [[ -n "${HOST_CDS_FNA}" ]]; then
    if [[ ! -f "${HOST_CDS_FNA}" ]]; then
        echo "ERROR: no host CDS file was find: ${HOST_CDS_FNA}" >&2
        exit 1
    fi
    echo "==> find CDS dile：${HOST_CDS_FNA}"
else
    echo "==> No CDS file was found, will use bedtools to get seq from host genome"
fi
if [[ ! -f "${PHAGE_FNA}" ]]; then
    echo "ERROR: NO phage genome file was found: ${PHAGE_FNA}" >&2
    exit 1
fi

#############################################
# check required program
#############################################

for cmd in makeblastdb blastn samtools bedtools seqkit; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: no command $cmd ,install module load first" >&2
        exit 1
    fi
done

#############################################
# setting output filename
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
# first: build blast DB
#############################################

echo "==> building BLAST database：${HOST_FNA}"
makeblastdb -in "${HOST_FNA}" -dbtype nucl -out "${BLAST_DB}" >/dev/null

#############################################
# second:blast phage genome to host genome
#############################################

echo "==> BLAST phage (${PHAGE_FNA}) to host (${HOST_FNA})"
blastn \
    -query "${PHAGE_FNA}" \
    -db "${BLAST_DB}" \
    -out "${BLAST_OUT}" \
    -outfmt 6 \
    -max_target_seqs 50 \
    -max_hsps 50

if [[ ! -s "${BLAST_OUT}" ]]; then
    echo "ERROR: no hit on blast, can't find phage region." >&2
    exit 1
fi

#############################################
# Third: from blast results to get prophage region -> PHAGE_BED
#
#   1) >3000 bp
#   2) calculate total_hit_len for each contig
#   3) pick the contig with highest total_hit_len
#   4) get min(start) and max(end) for that contig
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
# generate host clean regions bed
#############################################

echo "==> build genome index and get non-phage regions"

samtools faidx "${HOST_FNA}"
cut -f1,2 "${HOST_FNA}.fai" > "${HOST_GENOME_SIZES}"

bedtools complement -i "${PHAGE_BED}" -g "${HOST_GENOME_SIZES}" > "${HOST_CLEAN_REGIONS_BED}"

#############################################
# step5: get host genome without phage regions -> OUT_FNA_GENOMIC
# use bedtools getfasta, and merge all sequences to single header and single line sequence
#############################################

echo "==> generate ${OUT_FNA_GENOMIC}（host genome without phage region）"

bedtools getfasta \
    -fi "${HOST_FNA}" \
    -bed "${HOST_CLEAN_REGIONS_BED}" \
    -fo "${OUT_FNA_GENOMIC}.tmp"

# merge all sequences to single header and single line sequence
echo "==> merge all sequences to single header and single line sequence"
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
# step6: get all CDS -> BED
#############################################

echo "==> 從 GFF 轉成 CDS BED"

# check if is Bacteroides format（source column "Protein Homology"）
BACTEROIDES_FORMAT=0
head -1000 "${HOST_GFF}" > /tmp/gff_check.txt || true
if grep "Protein.*Homology.*CDS" /tmp/gff_check.txt > /dev/null 2>&1; then
    BACTEROIDES_FORMAT=1
fi
rm -f /tmp/gff_check.txt

if [[ ${BACTEROIDES_FORMAT} -eq 1 ]]; then
    echo "    detect Bacteroides GFF format（source is 'Protein Homology'）"
    
    # use special Bacteroides GFF to BED conversion
    # Bacteroides GFF format：Col1=seqid, Col2=Protein, Col3=Homology, Col4=CDS, Col5=start, Col6=end, Col7=., Col8=strand, Col9=phase, Col10=attributes
    awk '
    $2 == "Protein" && $3 == "Homology" && $4 == "CDS" {
        start = $5 - 1  # GFF is 1-based, BED need 0-based
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
        # GFF: col4, col5 is 1-based，BED need 0-based, 所以下限要 -1
        start = $4 - 1
        if (start < 0) start = 0
        printf("%s\t%d\t%s\t%s\t0\t%s\n", $1, start, $5, $9, $7)
    }
    ' "${HOST_GFF}" > "${HOST_CDS_BED}"
fi

#############################################
# 7.remove overlapped region of phage to CDS
#############################################

echo "==> remove CDS located in phage region p

bedtools subtract \
    -a "${HOST_CDS_BED}" \
    -b "${PHAGE_BED}" \
    > "${HOST_CDS_CLEAN_BED}"

#############################################
# 8. CDS -> Host_clean_CDS.fna
#############################################

if [[ -n "${HOST_CDS_FNA}" ]]; then
    echo "==> select phage region from CDS file -> ${OUT_FNA_CDS}"
    
    # extract CDS ID from HOST_CDS_CLEAN_BED (get rid of cds- for seqkit grep)
    grep -o "ID=[^;]*" "${HOST_CDS_CLEAN_BED}" | sed "s/ID=//" | sed "s/^cds-//" | sort -u > "${HOST_CLEAN_IDS}"
    
    # use awk to extract sequences from HOST_CDS_FNA based on IDs in HOST_CLEAN_IDS
    awk '
    BEGIN {
        # Read IDs into an array
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
    # if no CDS file, use bedtools to get sequences from genome
    echo "==> extract clean CDS seq from host genome -> ${OUT_FNA_CDS}"
    
    bedtools getfasta \
        -fi "${HOST_FNA}" \
        -bed "${HOST_CDS_CLEAN_BED}" \
        -s \
        -name \
        > "${OUT_FNA_CDS}"
fi

#############################################
# 8.： Host_clean.fnn
#############################################

echo "==> generate ${OUT_FNN}（only including host CDS nucleotide）"

bedtools getfasta \
    -fi "${HOST_FNA}" \
    -bed "${HOST_CDS_CLEAN_BED}" \
    -s \
    -name \
    > "${OUT_FNN}"

#############################################
# Step 9：according to CDS ID, get clean protein -> Host_clean.faa
#############################################

echo "==> generate ${OUT_FAA}（only host protein）"

# extract ID=xxx from 4 column of BED file, get rid of ID= and cds-
# if file exists, skip
if [[ ! -f "${HOST_CLEAN_IDS}" ]]; then
    grep -o "ID=[^;]*" "${HOST_CDS_CLEAN_BED}" | sed "s/ID=//" | sed "s/^cds-//" | sort -u > "${HOST_CLEAN_IDS}"
fi

# according to IDs to grep sequences from HOST_FAA
seqkit grep -f "${HOST_CLEAN_IDS}" "${HOST_FAA}" > "${OUT_FAA}"

echo "==> complete! output file："
echo "    ${OUT_FNA_GENOMIC}"
echo "    ${OUT_FNA_CDS}"
echo "    ${OUT_FNN}"
echo "    ${OUT_FAA}"
