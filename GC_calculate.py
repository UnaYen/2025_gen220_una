#!/usr/bin/env python3
"""
GC Content Calculator for FASTA sequences
Calculates GC percentage for each sequence and provides distribution analysis
"""

import sys
import os
import glob
from pathlib import Path
from collections import defaultdict

def calculate_gc_percentage(sequence):
    """
    Calculate GC percentage for a DNA sequence
    GC% = (G + C) / Total_Length * 100
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_length = len(sequence)
    
    if total_length == 0:
        return 0.0
    
    return (gc_count / total_length) * 100

def get_gc_bin(gc_percentage):
    """
    Determine which bin the GC percentage falls into (10%, 20%, 30%, etc.)
    Returns bin range as string (e.g., "10-20%", "20-30%")
    """
    bin_size = 10
    lower_bound = int(gc_percentage // bin_size) * bin_size
    upper_bound = lower_bound + bin_size
    return f"{lower_bound}-{upper_bound}%"

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
                    # Process previous sequence
                    if header is not None:
                        yield header, ''.join(sequence)
                    
                    # Start new sequence
                    header = line[1:]  # Remove '>' character
                    sequence = []
                else:
                    # Add to current sequence
                    if header is not None:
                        sequence.append(line)
            
            # Don't forget the last sequence
            if header is not None:
                yield header, ''.join(sequence)
    
    except IOError as e:
        print(f"❌ Error reading file {fasta_file}: {e}", file=sys.stderr)
        return

def analyze_gc_content(fasta_file):
    """
    Analyze GC content of all sequences in a FASTA file
    Returns:
        gc_distribution: dict with bin -> count
        individual_gc: list of (header, gc_percentage)
        overall_gc: overall GC percentage
    """
    gc_distribution = defaultdict(int)
    individual_gc = []
    total_gc_count = 0
    total_length = 0
    
    for header, sequence in parse_fasta(fasta_file):
        if not sequence:
            continue
        
        gc_percentage = calculate_gc_percentage(sequence)
        gc_bin = get_gc_bin(gc_percentage)
        
        gc_distribution[gc_bin] += 1
        individual_gc.append((header, gc_percentage, len(sequence)))
        
        # Calculate overall GC content
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        total_gc_count += gc_count
        total_length += len(sequence)
    
    # Calculate overall GC percentage
    overall_gc = (total_gc_count / total_length * 100) if total_length > 0 else 0.0
    
    return gc_distribution, individual_gc, overall_gc, total_length

def find_fasta_files(directory, pattern="*_CDS.fna"):
    """
    Find all FASTA files matching the pattern in the directory
    """
    file_path = os.path.join(directory, pattern)
    return sorted(glob.glob(file_path))

def print_report(filename, gc_distribution, individual_gc, overall_gc, total_length):
    """
    Print formatted report of GC content analysis
    """
    print(f"\n{'='*70}")
    print(f"GC Content Analysis Report")
    print(f"{'='*70}")
    print(f"File: {filename}")
    print(f"Total sequences: {len(individual_gc)}")
    print(f"Total length: {total_length:,} bp")
    print(f"Overall GC content: {overall_gc:.2f}%")
    print(f"\n{'-'*70}")
    
    # Print distribution histogram
    print(f"{'GC Content Bin':<20} {'Count':<10} {'Percentage':<15}")
    print(f"{'-'*70}")
    
    # Sort bins for display
    bin_order = sorted(gc_distribution.keys(), key=lambda x: int(x.split('-')[0]))
    
    total_sequences = len(individual_gc)
    for bin_range in bin_order:
        count = gc_distribution[bin_range]
        percentage = (count / total_sequences * 100) if total_sequences > 0 else 0
        bar = '█' * int(percentage / 2)  # Visual bar (2% per block)
        print(f"{bin_range:<20} {count:<10} {percentage:>6.2f}%  {bar}")
    
    print(f"{'-'*70}\n")

def main():
    """
    Main function to handle command-line arguments and file processing
    """
    if len(sys.argv) < 2:
        print("Usage: python3 GC_calculate.py [directory] [file_pattern]")
        print("\nExamples:")
        print("  python3 GC_calculate.py                    # Process core_gene directory")
        print("  python3 GC_calculate.py core_gene          # Specify directory")
        print("  python3 GC_calculate.py core_gene '*_CDS.fna'  # Custom pattern")
        print("  python3 GC_calculate.py /path/to/file.fna  # Single file")
        sys.exit(1)
    
    # Determine input
    input_arg = sys.argv[1]
    pattern = sys.argv[2] if len(sys.argv) > 2 else "*_CDS.fna"
    
    # Check if input is a file or directory
    if os.path.isfile(input_arg):
        # Single file input
        files_to_process = [input_arg]
    elif os.path.isdir(input_arg):
        # Directory input
        files_to_process = find_fasta_files(input_arg, pattern)
        if not files_to_process:
            print(f"❌ No files matching pattern '{pattern}' found in '{input_arg}'")
            sys.exit(1)
    else:
        print(f"❌ Invalid input: '{input_arg}' is neither a file nor a directory")
        sys.exit(1)
    
    # Process each file
    for fasta_file in files_to_process:
        print(f"📖 Processing: {fasta_file}")
        
        try:
            gc_distribution, individual_gc, overall_gc, total_length = analyze_gc_content(fasta_file)
            
            # Print report
            print_report(os.path.basename(fasta_file), gc_distribution, individual_gc, overall_gc, total_length)
            
            # Optional: Show top sequences by GC content (uncomment if needed)
            # print(f"Top 10 highest GC content sequences:")
            # sorted_gc = sorted(individual_gc, key=lambda x: x[1], reverse=True)
            # for i, (header, gc_pct, length) in enumerate(sorted_gc[:10], 1):
            #     print(f"  {i:2d}. {header[:50]:<50} {gc_pct:>6.2f}%  ({length} bp)")
            
        except Exception as e:
            print(f"❌ Error processing {fasta_file}: {e}", file=sys.stderr)
            continue

if __name__ == "__main__":
    main()
