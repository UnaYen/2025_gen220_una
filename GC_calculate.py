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

def print_report(filename, gc_distribution, individual_gc, overall_gc, total_length, detail_file=None):
    """
    Print formatted report of GC content analysis and optionally save to file
    """
    report_lines = []
    
    report_lines.append(f"\n{'='*70}")
    report_lines.append(f"GC Content Analysis Report")
    report_lines.append(f"{'='*70}")
    report_lines.append(f"File: {filename}")
    report_lines.append(f"Total sequences: {len(individual_gc)}")
    report_lines.append(f"Total length: {total_length:,} bp")
    report_lines.append(f"Overall GC content: {overall_gc:.2f}%")
    report_lines.append(f"\n{'-'*70}")
    
    # Print distribution histogram
    report_lines.append(f"{'GC Content Bin':<20} {'Count':<10} {'Percentage':<15}")
    report_lines.append(f"{'-'*70}")
    
    # Sort bins for display
    bin_order = sorted(gc_distribution.keys(), key=lambda x: int(x.split('-')[0]))
    
    total_sequences = len(individual_gc)
    for bin_range in bin_order:
        count = gc_distribution[bin_range]
        percentage = (count / total_sequences * 100) if total_sequences > 0 else 0
        bar = '█' * int(percentage / 2)  # Visual bar (2% per block)
        report_lines.append(f"{bin_range:<20} {count:<10} {percentage:>6.2f}%  {bar}")
    
    report_lines.append(f"{'-'*70}\n")
    
    # Print to console
    for line in report_lines:
        print(line)
    
    # Save to detail file if provided
    if detail_file:
        with open(detail_file, 'a') as f:
            for line in report_lines:
                f.write(line + '\n')
    
    return overall_gc

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
        print("\nOutput files:")
        print("  - GC_calculate_output_detail.txt  (detailed analysis)")
        print("  - GC_calculate_output_summary.txt (summary table)")
        sys.exit(1)
    
    # Output files
    detail_output = "GC_calculate_output_detail.txt"
    summary_output = "GC_calculate_output_summary.txt"
    
    # Clear output files if they exist
    open(detail_output, 'w').close()
    open(summary_output, 'w').close()
    
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
    
    # Store summary data
    summary_data = []
    
    # Process each file
    for fasta_file in files_to_process:
        print(f"📖 Processing: {fasta_file}")
        
        try:
            gc_distribution, individual_gc, overall_gc, total_length = analyze_gc_content(fasta_file)
            
            # Print report and save to detail file
            print_report(os.path.basename(fasta_file), gc_distribution, individual_gc, overall_gc, total_length, detail_output)
            
            # Store summary data
            summary_data.append((os.path.basename(fasta_file), overall_gc))
            
        except Exception as e:
            print(f"❌ Error processing {fasta_file}: {e}", file=sys.stderr)
            continue
    
    # Write summary table
    print(f"\n✓ Saving summary to {summary_output}...")
    with open(summary_output, 'w') as f:
        f.write(f"{'Filename':<50} {'Final GC%':<12}\n")
        f.write(f"{'-'*62}\n")
        for filename, gc_percent in sorted(summary_data):
            f.write(f"{filename:<50} {gc_percent:>10.2f}%\n")
    
    print(f"✓ Saving details to {detail_output}...")
    print(f"\n✅ Complete! Output files created:")
    print(f"   - {detail_output}")
    print(f"   - {summary_output}")

if __name__ == "__main__":
    main()
