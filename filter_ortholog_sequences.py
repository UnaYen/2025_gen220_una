#!/usr/bin/env python3
"""
according to Orthologs_all.txt, seperate core genes and other genes
Usage: python3 filter_ortholog_sequences.py <species_name> <input_dir> <output_dir>
"""

import sys
import os
from pathlib import Path

def load_ortholog_accessions(ortholog_file):
    """
    Read Orthologs_all.txt and return accession set of each species
    """
    ortholog_dict = {}
    current_species = None
    
    with open(ortholog_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # check if this is a species line
            if line.endswith(':'):
                current_species = line.rstrip(':')
                ortholog_dict[current_species] = set()
            elif current_species:
                # this is an accession
                ortholog_dict[current_species].add(line)
    
    return ortholog_dict

def extract_accession(header_line):
    """
    extract accession ID from FASTA header line
    deal with two formats:
    1. >ABR37731.1 ... (easy format)
    2. >lcl|JH724132.1_cds_EIY40997.1_1 [protein_id=EIY40997.1] ... (conplex format)
    """
    header_line = header_line.lstrip('>')
    
    # try extract from [protein_id=...]
    import re
    match = re.search(r'\[protein_id=([^\]]+)\]', header_line)
    if match:
        return match.group(1)
    
    # otherwise, take first word as accession
    accession = header_line.split()[0]
    return accession

def process_fasta_file(input_fasta, ortholog_accessions, output_core, output_other):
    """
    deal with a FASTA file and separate sequences into core and other based on ortholog_accessions
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
                    # write the previous sequence if exists
                    if current_header and current_accession:
                        sequence = ''.join(current_seq)
                        if current_accession in ortholog_accessions:
                            core_file.write(f"{current_header}\n{sequence}\n")
                            core_count += 1
                        else:
                            other_file.write(f"{current_header}\n{sequence}\n")
                            other_count += 1
                    
                    # get new header and accession
                    current_header = line
                    current_accession = extract_accession(line)
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # write the last sequence
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
    
    # make directory
    os.makedirs(output_dir, exist_ok=True)
    
    # access ortholog accessions
    ortholog_file = "Orthologs_all.txt"
    if not os.path.exists(ortholog_file):
        print(f"❌ can't find {ortholog_file}")
        sys.exit(1)
    
    print(f"📖 load {ortholog_file}...")
    ortholog_dict = load_ortholog_accessions(ortholog_file)
    
    if species_name not in ortholog_dict:
        print(f"❌ '{species_name}' not in Orthologs_all.txt")
        print(f"available species: {', '.join(sorted(ortholog_dict.keys()))}")
        sys.exit(1)
    
    ortholog_accessions = ortholog_dict[species_name]
    print(f"✓ load {len(ortholog_accessions)} {species_name} ortholog accessions")
    
    # find .faa and .fna files
    species_pattern = species_name.replace(':', '')
    
    potential_files = []
    
    # method 1: full species name
    potential_files.append(f"{input_dir}/{species_name}_protein.faa")
    potential_files.append(f"{input_dir}/{species_name}_cds_from_genomic.fna")
    
    # method 2: short name
    short_name = species_name.split('_')[-1]  # take last part as short name
    potential_files.append(f"{input_dir}/*{short_name}*protein.faa")
    potential_files.append(f"{input_dir}/*{short_name}*cds_from_genomic.fna")
    
    # check files
    faa_file = None
    fna_file = None
    
    for file_pattern in potential_files:
        if '*' in file_pattern:
            # use glob to find matching files
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
        print(f"\n❌ no {species_name} files found!")
        print(f"  search loaction: {input_dir}/")
        print(f"  expected .faa file: {faa_file if faa_file else '(未找到)'}")
        print(f"  expected .fna file: {fna_file if fna_file else '(未找到)'}")
        sys.exit(1)
    
    print(f"\n✓ file found:")
    print(f"  FAA: {faa_file}")
    print(f"  FNA: {fna_file}")
    
    # create output file names
    output_prefix = os.path.join(output_dir, species_name)
    core_faa = f"{output_prefix}_core_protein.faa"
    other_faa = f"{output_prefix}_other_protein.faa"
    core_fna = f"{output_prefix}_core_CDS.fna"
    other_fna = f"{output_prefix}_other_CDS.fna"
    
    # deal with .faa file
    print(f"\n🔄 deal with .faa file...")
    core_count_faa, other_count_faa = process_fasta_file(faa_file, ortholog_accessions, core_faa, other_faa)
    print(f"✓ Core proteins: {core_count_faa}")
    print(f"✓ Other proteins: {other_count_faa}")
    
    # deal with .fna file
    print(f"\n🔄 deal with .fna file...")
    core_count_fna, other_count_fna = process_fasta_file(fna_file, ortholog_accessions, core_fna, other_fna)
    print(f"✓ Core CDS: {core_count_fna}")
    print(f"✓ Other CDS: {other_count_fna}")
    
    # output summary
    print(f"\n{'='*60}")
    print(f"export file:")
    print(f"  {core_faa}")
    print(f"  {other_faa}")
    print(f"  {core_fna}")
    print(f"  {other_fna}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
