"""
Script to identify orthologs where Phocaeicola_vulgatus_ATCC has exactly one ortholog
in each of the other 9 bacterial species.
"""

import sys

# Define the target and other species
TARGET_SPECIES = "Phocaeicola_vulgatus_ATCC"
OTHER_SPECIES = [
    "Phocaeicola_dorei_HM719",
    "Bacteroides_thetaiotaomicron_VPI-5482",
    "Bacteroides_thetaiotaomicron_jmh43",
    "Bacteroides_intestinalis_APC919",
    "Bacteroides_uniformis_ATCC",
    "Alistipes_shahii_DA73",
    "Parabacteroides_distasonis_DA491",
    "Escherichia_coli_MG1655",
    "Escherichia_coli_BL21"
    ]

# Get input file from command line argument
if len(sys.argv) < 2:
    print("Usage: python find_ortholog.py <input_file>")
    print("Error: Input file is required")
    sys.exit(1)

input_file = sys.argv[1]

# Generate output file names based on input file
if input_file.endswith(".txt"):
    base_name = input_file[:-4]
else:
    base_name = input_file

output_file_detail = base_name + "_output_detail.txt"
output_file_summary = base_name + "_output_summary.txt"

# Dictionary to store results: {species: [accession_numbers]}
results = {species: [] for species in OTHER_SPECIES}

# Open output files for writing
with open(output_file_detail, "w") as detail_f, open(output_file_summary, "w") as summary_f:
    # Read and process the input file
    with open(input_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            # Split the tab-delimited line
            entries = line.strip().split("\t")
            
            # Count occurrences of target species
            target_count = sum(1 for entry in entries if entry.startswith(TARGET_SPECIES))
            
            # If there's exactly one Phocaeicola_vulgatus_ATCC
            if target_count == 1:
                # Check each other species
                for other_species in OTHER_SPECIES:
                    # Find entries from this species
                    species_entries = [entry for entry in entries if entry.startswith(other_species)]
                    
                    # If exactly one entry from this species exists
                    if len(species_entries) == 1:
                        # Extract the accession number (after the |)
                        accession = species_entries[0].split("|")[1]
                        results[other_species].append(accession)
                        msg = f"Line {line_num}: {other_species} -> {accession}"
                        detail_f.write(msg + "\n")
    
    # Write summary of results to summary file
    for species in OTHER_SPECIES:
        if results[species]:  # species has orthologs found
            summary_f.write(f"{species}:\n")
            for accession in results[species]:
                summary_f.write(f"{accession}\n")
            summary_f.write("\n")
