#!/bin/bash

# Example usage:
# INDIR=/data/rc003/atma/horizontalTransfer/results/L1/allSpeciesTogether FILE=allSpecies_L1.fasta MINCODONS=200 sbatch findOrfs.q

#SBATCH -A robinson
#SBATCH --qos=long
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=7-00:00
#SBATCH --mem=32GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au 

# Import BEDTools
module load BEDTools/2.25.0-foss-2015b

# Go to dir
cd $INDIR

# Use usearch to find all possible orfs 
# Output orfs in both nucleotide and amino acid form
usearch -fastx_findorfs "$FILE" -ntout "$FILE"_ORFnt.fasta -aaout "$FILE"_ORFaa.fasta -orfstyle 7 -mincodons "$MINCODONS"

# Scan candidates against known protein domains 
# using the Pfam-A.hmm database
hmmscan --domtblout "$FILE"_ORFaa.fasta.out /data/rc003/atma/horizontalTransfer/references/pfam/Pfam-A.hmm "$FILE"_ORFaa.fasta > "$FILE"_ORFaa.fasta.log

# pull out a BED file containing the location of the RT domains
# check that the RT domain validates as RVT_1 or RVT_3
# and that it is full-size (e.g. expect it to be around 255 amino acids,
# so set a restriction of >200 amino acids
# note: $20 and $21 are the envelope coordinates on the query protein seq
# so even if the alignment doesn't stretch as far, or scores low in some regions
# this 'envelope' spans the entire suspected RT domain
cat "$FILE"_ORFaa.fasta.out \
| awk '{if (($1=="RVT_1" || $1=="RVT_3") && ($21-$20>200)) print $4 "\t"  $20 "\t" $21}' \
> "$FILE"_RTdomain.bed

# retrieve the fasta using this BED file
# no need to use strand since they are all facing the same way
bedtools getfasta -fi "$FILE"_ORFaa.fasta -bed "$FILE"_RTdomain.bed -fo "$FILE"_RTdomain.fasta

