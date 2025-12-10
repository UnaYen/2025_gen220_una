"""
Script to identify orthologs where Phocaeicola_vulgatus_ATCC has exactly one ortholog
in each of the other 7 bacterial species.
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

# at least 4 of the following 6 species must have orthologs
CORE_SPECIES = [
    "Phocaeicola_dorei_HM719",
    "Bacteroides_thetaiotaomicron_VPI-5482",
    "Bacteroides_thetaiotaomicron_jmh43",
    "Bacteroides_intestinalis_APC919",
    "Bacteroides_uniformis_ATCC",
    "Parabacteroides_distasonis_DA491"
]

# 從命令行參數獲取輸入檔案名稱
if len(sys.argv) < 2:
    print("❌ usage: python3 find_ortholog_PV.py <input_file>")
    print("\nexample:")
    print("  python3 find_ortholog_PV.py Orthogroups_clean.txt")
    sys.exit(1)

input_file = sys.argv[1]
print(f"📖 Read input file: {input_file}")

# Dictionary to store results for each PV accession
# Format: {pv_accession: {"core_species": set(), "orthologs": {species: [accessions]}}}
pv_accession_data = {}

# Read and process the Orthogroups.txt file
with open(input_file, "r") as f:
    for line_num, line in enumerate(f, 1):
        # Split the tab-delimited line
        entries = line.strip().split("\t")
        
        # Count occurrences of target species
        target_count = sum(1 for entry in entries if entry.startswith(TARGET_SPECIES))
        
        # If there's exactly one Phocaeicola_vulgatus_ATCC
        if target_count == 1:
            # get Phocaeicola_vulgatus_ATCC accession number
            pv_entry = [entry for entry in entries if entry.startswith(TARGET_SPECIES)][0]
            pv_accession = pv_entry.split("|")[1]
            
            # initialize data structure for this PV accession if not exists
            if pv_accession not in pv_accession_data:
                pv_accession_data[pv_accession] = {
                    "core_species": set(),
                    "orthologs": {species: [] for species in CORE_SPECIES}
                }
            
            # check each core species
            for core_species in CORE_SPECIES:
                # Find entries from this species
                species_entries = [entry for entry in entries if entry.startswith(core_species)]
                
                # If exactly one entry from this species exists
                if len(species_entries) == 1:
                    accession = species_entries[0].split("|")[1]
                    pv_accession_data[pv_accession]["core_species"].add(core_species)
                    # Avoid duplicates
                    if not pv_accession_data[pv_accession]["orthologs"][core_species]:
                        pv_accession_data[pv_accession]["orthologs"][core_species].append(accession)
                        print(f"Line {line_num}: PV {pv_accession} + {core_species} -> {accession}")

# Print summary of results
print("\n" + "="*80)
print("SUMMARY OF RESULTS")
print("="*80)

# 篩選符合條件的 PV accession（至少 4 個核心物種）
pv_ortholog_dict = {}
for pv_accession, data in pv_accession_data.items():
    core_species_set = data["core_species"]
    core_count = len(core_species_set)
    
    if core_count >= 4:
        pv_ortholog_dict[pv_accession] = {
            "core_species_count": core_count,
            "orthologs": {sp: data["orthologs"][sp][0] if data["orthologs"][sp] else None 
                         for sp in core_species_set}
        }

print(f"\n找到 {len(pv_ortholog_dict)} 個符合條件的 Phocaeicola_vulgatus_ATCC 正交體")
print("(需要在核心物種中至少 4 個有正交體)\n")

for pv_accession in sorted(pv_ortholog_dict.keys()):
    data = pv_ortholog_dict[pv_accession]
    core_count = data['core_species_count']
    orthologs = data['orthologs']
    
    print(f"PV Accession: {pv_accession}")
    print(f"  核心物種正交體: {core_count}/6")
    for species, accession in orthologs.items():
        print(f"    - {species}: {accession}")
    print()

# Save results to file
output_file = "PV_ortholog.txt"
with open(output_file, "w") as f:
    f.write("="*80 + "\n")
    f.write("PHOCAEICOLA_VULGATUS ORTHOLOGS\n")
    f.write("(符合: 核心物種中至少 4 個有正交體)\n")
    f.write("="*80 + "\n\n")
    
    f.write(f"總計: {len(pv_ortholog_dict)} 個\n\n")
    
    for pv_accession in sorted(pv_ortholog_dict.keys()):
        data = pv_ortholog_dict[pv_accession]
        f.write(f"PV Accession: {pv_accession}\n")
        f.write(f"  核心物種正交體: {data['core_species_count']}/6\n")
        f.write(f"  正交體:\n")
        for species, accession in data['orthologs'].items():
            f.write(f"    - {species}: {accession}\n")
        f.write("\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("完整字典:\n")
    f.write("="*80 + "\n")
    f.write(str(pv_ortholog_dict))

print(f"\n✓ 結果已保存至 {output_file}")

# Print the PV orthologs dictionary
print("\n" + "="*80)
print("PHOCAEICOLA_VULGATUS ORTHOLOGS DICTIONARY")
print("="*80)
print(pv_ortholog_dict)
