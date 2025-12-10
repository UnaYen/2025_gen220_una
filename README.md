--clean_host.sh--
find out phage region in host genome, and get rid of prophage region, then output
1) Host_clean_genomic.fna  (including header and merged sequence)
2) Host_clean_CDS.fna      (including host CDS nucleotide)
3) Host_clean.fnn          (including host CDS nucleotide)
4) Host_clean.faa          (including host protein)

usage:
#   bash clean_host.sh -input HOST_PREFIX PHAGE_FNA -output all
Parameters:
 "  -input       HOST_FNA and PHAGE_FNA files (required)"
 "  -minlen      Minimum alignment length in bp (default: 3000)"
 "  -output      Output mode (default: all)"