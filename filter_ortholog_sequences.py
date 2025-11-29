#!/usr/bin/env python3
"""
根據 Orthologs_all.txt 將序列分類為 core genes 和 other genes
Usage: python3 filter_ortholog_sequences.py <species_name> <input_dir> <output_dir>
"""

import sys
import os
from pathlib import Path

def load_ortholog_accessions(ortholog_file):
    """
    讀取 Orthologs_all.txt 並返回每個物種的 accession set
    """
    ortholog_dict = {}
    current_species = None
    
    with open(ortholog_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # 檢查是否是物種名稱行（以冒號結尾）
            if line.endswith(':'):
                current_species = line.rstrip(':')
                ortholog_dict[current_species] = set()
            elif current_species:
                # 這是一個 accession
                ortholog_dict[current_species].add(line)
    
    return ortholog_dict

def extract_accession(header_line):
    """
    從 FASTA header 提取 accession ID
    處理兩種格式:
    1. >ABR37731.1 ... (簡單格式)
    2. >lcl|JH724132.1_cds_EIY40997.1_1 [protein_id=EIY40997.1] ... (複雜格式)
    """
    header_line = header_line.lstrip('>')
    
    # 首先嘗試從 [protein_id=...] 提取
    import re
    match = re.search(r'\[protein_id=([^\]]+)\]', header_line)
    if match:
        return match.group(1)
    
    # 否則取第一個 token
    accession = header_line.split()[0]
    return accession

def process_fasta_file(input_fasta, ortholog_accessions, output_core, output_other):
    """
    處理單個 FASTA 檔案，根據 accession 分類
    """
    core_count = 0
    other_count = 0
    
    with open(input_fasta, 'r') as infile:
        with open(output_core, 'w') as core_file, open(output_other, 'w') as other_file:
            current_header = None
            current_accession = None
            current_seq = []
            
            for line in infile:
                line = line.rstrip('\n')
                
                if line.startswith('>'):
                    # 寫入上一個序列
                    if current_header and current_accession:
                        sequence = ''.join(current_seq)
                        if current_accession in ortholog_accessions:
                            core_file.write(f"{current_header}\n{sequence}\n")
                            core_count += 1
                        else:
                            other_file.write(f"{current_header}\n{sequence}\n")
                            other_count += 1
                    
                    # 提取新的 header 和 accession
                    current_header = line
                    current_accession = extract_accession(line)
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # 寫入最後一個序列
            if current_header and current_accession:
                sequence = ''.join(current_seq)
                if current_accession in ortholog_accessions:
                    core_file.write(f"{current_header}\n{sequence}\n")
                    core_count += 1
                else:
                    other_file.write(f"{current_header}\n{sequence}\n")
                    other_count += 1
    
    return core_count, other_count

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 filter_ortholog_sequences.py <species_name> [input_dir] [output_dir]")
        print("\nExample: python3 filter_ortholog_sequences.py Phocaeicola_dorei_HM719")
        print("         python3 filter_ortholog_sequences.py 'Bacteroides_thetaiotaomicron_VPI-5482'")
        sys.exit(1)
    
    species_name = sys.argv[1]
    input_dir = sys.argv[2] if len(sys.argv) > 2 else "genome_seq"
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "core_gene"
    
    # 建立輸出目錄
    os.makedirs(output_dir, exist_ok=True)
    
    # 載入 ortholog accessions
    ortholog_file = "Orthologs_all.txt"
    if not os.path.exists(ortholog_file):
        print(f"❌ 找不到 {ortholog_file}")
        sys.exit(1)
    
    print(f"📖 載入 {ortholog_file}...")
    ortholog_dict = load_ortholog_accessions(ortholog_file)
    
    if species_name not in ortholog_dict:
        print(f"❌ 物種 '{species_name}' 不在 Orthologs_all.txt 中")
        print(f"可用的物種: {', '.join(sorted(ortholog_dict.keys()))}")
        sys.exit(1)
    
    ortholog_accessions = ortholog_dict[species_name]
    print(f"✓ 已載入 {len(ortholog_accessions)} 個 {species_name} ortholog accessions")
    
    # 找到相應的 .faa 和 .fna 檔案
    # 使用更靈活的方式尋找檔案
    species_pattern = species_name.replace(':', '')
    
    # 建立可能的檔案名稱列表
    potential_files = []
    
    # 方式 1: 完整物種名稱
    potential_files.append(f"{input_dir}/{species_name}_protein.faa")
    potential_files.append(f"{input_dir}/{species_name}_cds_from_genomic.fna")
    
    # 方式 2: 以底線分隔的物種名稱
    short_name = species_name.split('_')[-1]  # 取最後一個部分
    potential_files.append(f"{input_dir}/*{short_name}*protein.faa")
    potential_files.append(f"{input_dir}/*{short_name}*cds_from_genomic.fna")
    
    # 實際檢查檔案
    faa_file = None
    fna_file = None
    
    for file_pattern in potential_files:
        if '*' in file_pattern:
            # 使用 glob 尋找
            import glob
            matches = glob.glob(file_pattern)
            if matches:
                if 'protein' in file_pattern and not faa_file:
                    faa_file = matches[0]
                elif 'cds' in file_pattern and not fna_file:
                    fna_file = matches[0]
        else:
            if 'protein' in file_pattern and os.path.exists(file_pattern):
                faa_file = file_pattern
            elif 'cds' in file_pattern and os.path.exists(file_pattern):
                fna_file = file_pattern
    
    if not faa_file or not fna_file:
        print(f"\n❌ 找不到 {species_name} 的序列檔案")
        print(f"  尋找位置: {input_dir}/")
        print(f"  預期的 .faa 檔案: {faa_file if faa_file else '(未找到)'}")
        print(f"  預期的 .fna 檔案: {fna_file if fna_file else '(未找到)'}")
        sys.exit(1)
    
    print(f"\n✓ 找到檔案:")
    print(f"  FAA: {faa_file}")
    print(f"  FNA: {fna_file}")
    
    # 建立輸出檔案名稱
    output_prefix = os.path.join(output_dir, species_name)
    core_faa = f"{output_prefix}_core_protein.faa"
    other_faa = f"{output_prefix}_other_protein.faa"
    core_fna = f"{output_prefix}_core_CDS.fna"
    other_fna = f"{output_prefix}_other_CDS.fna"
    
    # 處理 .faa 檔案
    print(f"\n🔄 處理 .faa 檔案...")
    core_count_faa, other_count_faa = process_fasta_file(faa_file, ortholog_accessions, core_faa, other_faa)
    print(f"✓ Core proteins: {core_count_faa}")
    print(f"✓ Other proteins: {other_count_faa}")
    
    # 處理 .fna 檔案
    print(f"\n🔄 處理 .fna 檔案...")
    core_count_fna, other_count_fna = process_fasta_file(fna_file, ortholog_accessions, core_fna, other_fna)
    print(f"✓ Core CDS: {core_count_fna}")
    print(f"✓ Other CDS: {other_count_fna}")
    
    # 統計輸出
    print(f"\n{'='*60}")
    print(f"輸出檔案:")
    print(f"  {core_faa}")
    print(f"  {other_faa}")
    print(f"  {core_fna}")
    print(f"  {other_fna}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
