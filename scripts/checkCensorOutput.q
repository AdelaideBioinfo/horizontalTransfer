#!/bin/bash

# Want to keep only the nucleotide sequences that CENSOR recognises as L1s or BovBs
# And not some other repeat (e.g. CR1 or hAT)
# To check L1s, set:
# ELEMENT=L1 
# QUERY=/data/rc003/atma/horizontalTransfer/queries/1826_L1_and_Tx1_repbase_headers.txt
# To check BovBs, set:
# ELEMENT=BovB 
# QUERY=/data/rc003/atma/horizontalTransfer/queries/37_BovB_headers.txt

# Example usage:
# SPECIES=Yarrowia.lipolytica FILE=Yarrowia.lipolytica_L1_combined.fasta GENOME=/data/rc003/atma/horizontalTransfer/genomes/fungi/Yarrowia.lipolytica/Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna ELEMENT=L1 QUERY=/data/rc003/atma/horizontalTransfer/queries/1826_L1_and_Tx1_repbase_headers.txt sbatch checkCensorOutput.q

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --time=0:10:00 
#SBATCH --mem=1GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au      

# Load the necessary modules
module load BEDTools/2.25.0-foss-2015b

# Go to CENSOR results directory
cd /data/rc003/atma/horizontalTransfer/results/$ELEMENT/censored/

# First, extract map file lines in the right orientation (d, not c)
# Print the header (which contains coords and strand) and column 4 (hit name)
# Put into BED-like format
cat "$FILE".map \
| awk '{if ($7=="d") print $0 }' \
| awk '{print $1 "\t" $4}' \
| sed 's/:/\t/g' \
| sed 's/(-)/.rev/g' \
| sed 's/(+)/.fwd/g' \
| awk '{gsub(/-/,"\t",$2); print}' \
| sed 's/.rev/\t-/g' \
| sed 's/.fwd/\t+/g' \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" "1" "\t" $4}' \
> "$FILE".rearranged.tmp

# Sort the file by chr/scaffold then coordinates
bedtools sort -i "$FILE".rearranged.tmp > "$FILE".sorted.tmp

# Merge overlapping hits
bedtools merge -s -i "$FILE".sorted.tmp -c 4 -o collapse > "$FILE".merged.tmp

# Rearrange columns 
awk '{print $5 "\t" $1 "\t" $2 "\t" $3 "\t" $4}' "$FILE".merged.tmp > "$FILE".mergedbed.tmp

# Extract sequences which identified as L1s (or, BovBs)
# For L1s: download from Repbase all potential L1 sequences
# This includes everything in the L1 and Tx1 subgroups (1826 seqs total)
# Extract all the header names -> 1826_L1_and_Tx1_repbase_headers.txt
# Compare this to the hit names in the merged.txt file
# and then rearrange into BED format again
awk 'NR==FNR { a[$1] = $1; next} { for (k in a) if ($1 ~ a[k]) { print $0; break } }' "$QUERY" "$FILE".mergedbed.tmp \
| awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5}' \
> "$SPECIES"_"$ELEMENT"_censored.bed

# Extract FASTA from L1-only merged BED file
bedtools getfasta -s -fi "$GENOME" -bed "$SPECIES"_"$ELEMENT"_censored.bed -fo "$SPECIES"_"$ELEMENT"_censored.tmp

# Sort sequences by length
usearch -sortbylength "$SPECIES"_"$ELEMENT"_censored.tmp -fastaout "$SPECIES"_"$ELEMENT"_censored.fasta

# Remove unnecessary files
rm *.tmp
