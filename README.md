--download.txt--
-
>/main/genome_seq/download.txt

list of bacterial and phage genome downloaded list, except for phage PV01. PV01 genome data was in /genome_seq directory.

--clean_host.sh--
-
input one bacterial genomic.fna file and phage genomic.fna file, find out phage region in host genome, and get rid of prophage region, then output
1) Host_clean_genomic.fna  (including header and merged sequence)
2) Host_clean_CDS.fna      (including host CDS nucleotide)
3) Host_clean.fnn          (including host CDS nucleotide)
4) Host_clean.faa          (including host protein)

usage:

    bash clean_host.sh -input HOST_PREFIX PHAGE_FNA -output all

Parameters:

 -input    >HOST_FNA and PHAGE_FNA files (required)
 
 -minlen    >Minimum alignment length in bp (default: 3000)
 
 -output    >Output mode (default: all)

 --find_ortholog.py--
 -
>get orthologs list for bacteria except for P.vulgatus

usage:

    python ./find_ortholog.py Orthogroups_clean.txt

 1. input Ortholog_clean.txt, which got from Orthovenn3 output.
 2. orthologs difinition: if there are  only onc PV protein and other species's protien show up in same clade, consider that protein to be an ortholos to PV species.

output file: Orthogroups_clean_output_summary.txt

 --find_ortholog_PV.py--
-
>get PV core_gene list

usage:

python3 find_ortholog_PV.py Orthogroups_clean.txt

1. input Ortholog_clean.txt, which got from Orthovenn3 output.

2. PV_core gene difinition: 
set 6 related bacteria to be "core_species", if there are >=4 species have orthologs to same PV protein, consider that PV protein as PC_core gene/protien.

output file: PV_ortholog_summary.txt

---ortholog_all.txt--
-
>add two file content manually

Orthogroups_clean_output_summary.txt + PV_ortholog_summary.txt

 --find_ortholog_sequence---
-
>based on Ortholog list, filtered which gene belongs to core gene set and which belongs to accessary (other) gene set then output fasta file .

usage:

1. to batch process all data list in download.txt:

       bash batch_filter_orthologs.sh
2. or to process single file:

       python3 filter_ortholog_sequences.py "$species" "genome_seq" "core_gene"

output file:
1) Host_name_core_faa/fna
2) Host_name_other_faa/fna
main/core_gene

--GC_calculate.py--
-
>process for CDS_file, calculated each sequence and use 10% bin size, and also process final GC% for entire sequences.

usage:

    python3 GC_calculate.py [directory] [file_pattern]
      ex1: python3 GC_calculate.py core_gene '*_CDS.fna'
      ex2: python3 GC_calculate.py /path/to/file.fna
output file:
1. GC_calculate_output_detail.txt (with GC content bin distribution)
2. GC_calculate_output_summary.txt (final_GC% of each input file)

--codon_usage.py--
-
>Calculate codon usage and compare to other genome

usage:

    python3 codon_usage.py [directory] [compare_file]

output file:
1. count_hostname rscu_hostname file
   (main/codon_count directory)
2. codon_usage_cosine_similarity.txt
3. Instead of generate new odon_usage_cosine_similarity.txt, if the output file is already exist, generate new column in exiseted file for new output.
 
