#!/usr/bin/env python3
"""
Codon Usage Analysis Script
- Calculate codon frequency for each sequence in FASTA files
- Calculate RSCU (Relative Synonymous Codon Usage)
- Compare codon usage patterns using cosine similarity
"""

import sys
import os
import glob
from collections import defaultdict
from pathlib import Path
import math

# Standard Genetic Code Table 11 (Bacterial)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Create reverse mapping: amino acid -> list of codons
AA_TO_CODONS = defaultdict(list)
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS[aa].append(codon)

def parse_fasta(fasta_file):
    """
    Parse FASTA file and yield (header, sequence) tuples
    """
    header = None
    sequence = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                
                if line.startswith('>'):
                    if header is not None:
                        yield header, ''.join(sequence)
                    header = line[1:]
                    sequence = []
                else:
                    if header is not None:
                        sequence.append(line.upper())
            
            if header is not None:
                yield header, ''.join(sequence)
    
    except IOError as e:
        print(f"❌ Error reading file {fasta_file}: {e}", file=sys.stderr)
        return

def extract_codons(sequence):
    """
    Extract codons from sequence (3-base groups)
    Skip incomplete codons at the end
    """
    codons = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3 and all(base in 'ATCG' for base in codon):
            codons.append(codon)
    return codons

def count_codons_in_file(fasta_file):
    """
    Count all codons in a FASTA file
    Returns dict: codon -> count
    """
    codon_counts = defaultdict(int)
    total_sequences = 0
    total_codons = 0
    
    for header, sequence in parse_fasta(fasta_file):
        if not sequence:
            continue
        
        codons = extract_codons(sequence)
        for codon in codons:
            codon_counts[codon] += 1
            total_codons += 1
        
        total_sequences += 1
    
    return dict(codon_counts), total_sequences, total_codons

def calculate_rscu(codon_counts):
    """
    Calculate RSCU (Relative Synonymous Codon Usage)
    RSCU = (codon count) / (average count of synonymous codons)
    """
    rscu = {}
    
    # Group codons by amino acid
    for aa, codons in AA_TO_CODONS.items():
        if aa == '*':  # Skip stop codons
            continue
        
        # Get counts for all synonymous codons
        counts = [codon_counts.get(codon, 0) for codon in codons]
        total_count = sum(counts)
        
        if total_count == 0:
            # If no synonymous codons found, RSCU = 0
            for codon in codons:
                rscu[codon] = 0.0
        else:
            # Average count per synonymous codon
            avg_count = total_count / len(codons)
            
            # Calculate RSCU for each codon
            for codon in codons:
                codon_count = codon_counts.get(codon, 0)
                if avg_count > 0:
                    rscu[codon] = codon_count / avg_count
                else:
                    rscu[codon] = 0.0
    
    return rscu

def cosine_similarity(rscu1, rscu2):
    """
    Calculate cosine similarity between two RSCU vectors
    Cosine similarity = (dot product) / (magnitude1 * magnitude2)
    """
    # Get all unique codons
    all_codons = set(rscu1.keys()) | set(rscu2.keys())
    
    dot_product = 0.0
    mag1_squared = 0.0
    mag2_squared = 0.0
    
    for codon in all_codons:
        val1 = rscu1.get(codon, 0.0)
        val2 = rscu2.get(codon, 0.0)
        
        dot_product += val1 * val2
        mag1_squared += val1 * val1
        mag2_squared += val2 * val2
    
    mag1 = math.sqrt(mag1_squared)
    mag2 = math.sqrt(mag2_squared)
    
    if mag1 == 0 or mag2 == 0:
        return 0.0
    
    return dot_product / (mag1 * mag2)

def find_fasta_files(directory, pattern="*.fna"):
    """
    Find all FASTA files matching pattern in directory
    """
    file_path = os.path.join(directory, pattern)
    return sorted(glob.glob(file_path))

def save_codon_counts(filename, codon_counts, output_file):
    """
    Save codon counts to file
    Format: codon\tcount
    """
    with open(output_file, 'w') as f:
        f.write("Codon\tCount\n")
        # Sort by codon for consistency
        for codon in sorted(codon_counts.keys()):
            f.write(f"{codon}\t{codon_counts[codon]}\n")

def save_rscu(filename, rscu, output_file):
    """
    Save RSCU values to file
    Format: codon\trscu_value
    """
    with open(output_file, 'w') as f:
        f.write("Codon\tRSCU\tAminoAcid\n")
        # Sort by codon for consistency
        for codon in sorted(rscu.keys()):
            aa = CODON_TABLE.get(codon, '?')
            f.write(f"{codon}\t{rscu[codon]:.4f}\t{aa}\n")

def main():
    """
    Main function to analyze codon usage
    """
    if len(sys.argv) < 2:
        print("Usage: python3 codon_usage.py [directory] [compare_file]")
        print("\nExamples:")
        print("  python3 codon_usage.py core_gene")
        print("  python3 codon_usage.py core_gene 'Phocaeicola_vulgatus_ATCC_core_CDS.fna'")
        print("\nOutput:")
        print("  - codon_count/count_*.txt: Codon frequency for each file")
        print("  - codon_count/rscu_*.txt: RSCU values for each file")
        print("  - codon_usage_cosine_similarity.txt: Similarity comparison (if compare_file specified)")
        sys.exit(1)
    
    directory = sys.argv[1]
    compare_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Create output directory for codon counts
    codon_count_dir = "codon_count"
    os.makedirs(codon_count_dir, exist_ok=True)
    
    # Verify directory exists
    if not os.path.isdir(directory):
        print(f"❌ Directory not found: {directory}")
        sys.exit(1)
    
    # Find FASTA files
    fasta_files = find_fasta_files(directory, "*.fna")
    if not fasta_files:
        print(f"❌ No .fna files found in {directory}")
        sys.exit(1)
    
    # If reference file is provided and it's not in the search directory, add it to the list
    if compare_file and os.path.isfile(compare_file):
        if compare_file not in fasta_files:
            fasta_files.append(compare_file)
    
    print(f"Found {len(fasta_files)} FASTA files\n")
    
    # Process each file
    all_rscu = {}
    reference_rscu = None
    reference_name = None
    
    for fasta_file in fasta_files:
        basename = os.path.basename(fasta_file)
        print(f"📖 Processing: {basename}")
        
        try:
            # Count codons
            codon_counts, num_seqs, num_codons = count_codons_in_file(fasta_file)
            
            # Calculate RSCU
            rscu = calculate_rscu(codon_counts)
            all_rscu[basename] = rscu
            
            # Save results to codon_count directory
            count_output = os.path.join(codon_count_dir, f"count_{basename.replace('.fna', '.txt')}")
            rscu_output = os.path.join(codon_count_dir, f"rscu_{basename.replace('.fna', '.txt')}")
            
            save_codon_counts(basename, codon_counts, count_output)
            save_rscu(basename, rscu, rscu_output)
            
            print(f"  ✓ Sequences: {num_seqs}")
            print(f"  ✓ Total codons: {num_codons:,}")
            print(f"  ✓ Unique codons: {len(codon_counts)}")
            print(f"  ✓ Saved: {count_output}, {rscu_output}")
            
    # Track reference file for comparison
            if compare_file:
                # Check both full path and basename
                compare_basename = os.path.basename(compare_file)
                if basename == compare_file or basename == compare_basename:
                    reference_rscu = rscu
                    reference_name = basename
            
        except Exception as e:
            print(f"  ❌ Error: {e}", file=sys.stderr)
            continue
        
        print()
    
    # Calculate cosine similarity if reference file specified
    if compare_file:
        if reference_rscu is None:
            print(f"❌ Reference file '{compare_file}' not found")
            sys.exit(1)
        
        print(f"\n{'='*70}")
        print(f"Cosine Similarity Comparison (Reference: {reference_name})")
        print(f"{'='*70}\n")
        
        similarity_output = "codon_usage_cosine_similarity.txt"
        similarities = []
        
        with open(similarity_output, 'w') as f:
            # Write reference file information as first line
            f.write(f"Reference_File\t{reference_name}\n")
            f.write(f"Filename\tCosine_Similarity\n")
            
            for filename in sorted(all_rscu.keys()):
                if filename == reference_name:
                    continue
                
                rscu = all_rscu[filename]
                similarity = cosine_similarity(reference_rscu, rscu)
                similarities.append((filename, similarity))
                
                f.write(f"{filename}\t{similarity:.6f}\n")
                print(f"{filename:<50} {similarity:.6f}")
        
        print(f"\n✓ Saved: {similarity_output}")
        
        # Summary statistics
        if similarities:
            sims = [s[1] for s in similarities]
            avg_sim = sum(sims) / len(sims)
            max_sim = max(sims)
            min_sim = min(sims)
            
            print(f"\nSummary Statistics:")
            print(f"  Average Similarity: {avg_sim:.6f}")
            print(f"  Max Similarity: {max_sim:.6f}")
            print(f"  Min Similarity: {min_sim:.6f}")

if __name__ == "__main__":
    main()
