#!/usr/bin/env python3
"""
Aggregate codon count and RSCU files from codon_count directory into summary files.
Creates two summary files:
1. codon_count_summary.txt - All count files merged with filenames as columns
2. codon_rscu_summary.txt - All RSCU files merged with filenames as columns
"""

import os
import re
from collections import defaultdict

# Configuration
codon_dir = "codon_count"
output_dir = "."

def extract_sample_name(filename):
    """Extract sample name from count/rscu filename"""
    # Remove prefix (count_ or rscu_) and suffix (.txt)
    name = filename.replace("count_", "").replace("rscu_", "").replace(".txt", "")
    return name

def read_count_file(filepath):
    """Read count file and return dict {codon: count}"""
    data = {}
    try:
        with open(filepath, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    codon = parts[0]
                    count = parts[1]
                    data[codon] = count
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}")
    return data

def read_rscu_file(filepath):
    """Read RSCU file and return dict {codon: (rscu, amino_acid)}"""
    data = {}
    try:
        with open(filepath, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    codon = parts[0]
                    rscu = parts[1]
                    amino_acid = parts[2]
                    data[codon] = (rscu, amino_acid)
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}")
    return data

def get_all_codons():
    """Return sorted list of all valid codons"""
    bases = ['A', 'T', 'C', 'G']
    codons = []
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                codons.append(b1 + b2 + b3)
    return sorted(codons)

def aggregate_counts(codon_dir):
    """Aggregate all count files"""
    count_files = sorted([f for f in os.listdir(codon_dir) if f.startswith("count_")])
    
    if not count_files:
        print("No count files found!")
        return None, None
    
    print(f"Found {len(count_files)} count files")
    
    # Read all files
    all_counts = {}
    sample_names = []
    
    for count_file in count_files:
        filepath = os.path.join(codon_dir, count_file)
        sample_name = extract_sample_name(count_file)
        sample_names.append(sample_name)
        all_counts[sample_name] = read_count_file(filepath)
        print(f"  ✓ {sample_name}")
    
    return all_counts, sample_names

def aggregate_rscu(codon_dir):
    """Aggregate all RSCU files"""
    rscu_files = sorted([f for f in os.listdir(codon_dir) if f.startswith("rscu_")])
    
    if not rscu_files:
        print("No RSCU files found!")
        return None, None
    
    print(f"Found {len(rscu_files)} RSCU files")
    
    # Read all files
    all_rscu = {}
    sample_names = []
    
    for rscu_file in rscu_files:
        filepath = os.path.join(codon_dir, rscu_file)
        sample_name = extract_sample_name(rscu_file)
        sample_names.append(sample_name)
        all_rscu[sample_name] = read_rscu_file(filepath)
        print(f"  ✓ {sample_name}")
    
    return all_rscu, sample_names

def write_count_summary(all_counts, sample_names, output_dir):
    """Write aggregated count summary"""
    output_file = os.path.join(output_dir, "codon_count_summary.txt")
    codons = get_all_codons()
    
    print(f"\nWriting count summary to {output_file}...")
    
    with open(output_file, 'w') as f:
        # Write header
        header = ["Codon"] + sample_names
        f.write('\t'.join(header) + '\n')
        
        # Write data rows
        for codon in codons:
            row = [codon]
            for sample_name in sample_names:
                count = all_counts[sample_name].get(codon, "0")
                row.append(count)
            f.write('\t'.join(row) + '\n')
    
    print(f"✓ Count summary created: {output_file}")

def write_rscu_summary(all_rscu, sample_names, output_dir):
    """Write aggregated RSCU summary"""
    output_file = os.path.join(output_dir, "codon_rscu_summary.txt")
    codons = get_all_codons()
    
    print(f"\nWriting RSCU summary to {output_file}...")
    
    with open(output_file, 'w') as f:
        # Write header with codon info and all samples
        header = ["Codon", "AminoAcid"] + sample_names
        f.write('\t'.join(header) + '\n')
        
        # Write data rows
        for codon in codons:
            row = [codon]
            
            # Get amino acid from first file that has it
            amino_acid = ""
            for sample_name in sample_names:
                if codon in all_rscu[sample_name]:
                    amino_acid = all_rscu[sample_name][codon][1]
                    break
            row.append(amino_acid)
            
            # Add RSCU values for each sample
            for sample_name in sample_names:
                rscu = all_rscu[sample_name].get(codon, ("", ""))[0]
                row.append(rscu)
            
            f.write('\t'.join(row) + '\n')
    
    print(f"✓ RSCU summary created: {output_file}")

def main():
    print(f"{'='*70}")
    print("Codon Summary Aggregation")
    print(f"{'='*70}\n")
    
    # Aggregate count files
    print("Reading count files...")
    all_counts, count_samples = aggregate_counts(codon_dir)
    
    if all_counts:
        write_count_summary(all_counts, count_samples, output_dir)
    
    # Aggregate RSCU files
    print("\nReading RSCU files...")
    all_rscu, rscu_samples = aggregate_rscu(codon_dir)
    
    if all_rscu:
        write_rscu_summary(all_rscu, rscu_samples, output_dir)
    
    print(f"\n{'='*70}")
    print("✓ Aggregation complete!")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
