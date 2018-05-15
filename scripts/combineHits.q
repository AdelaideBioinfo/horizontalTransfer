#!/bin/bash

# Example usage:
# SPECIES=Yarrowia.lipolytica ELEMENT=L1 LASTZFILE=/data/rc003/atma/horizontalTransfer/results/L1/lastzed/Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna_Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna_L1_ORFp.fasta_nucl_seqs.fasta_lastz.bed TBLASTNFILE=/data/rc003/atma/horizontalTransfer/results/L1/tblastned/Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna_L1_ORFp.fasta_merged.bed GENOME=/data/rc003/atma/horizontalTransfer/genomes/fungi/Yarrowia.lipolytica/Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna RESULTSDIR=/data/rc003/atma/horizontalTransfer/results/L1/combined sbatch combineHits.q

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --time=0:10:00 
#SBATCH --mem=1GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au      

module load BEDTools/2.25.0-foss-2015b

cat $TBLASTNFILE $LASTZFILE > "$SPECIES"_"$ELEMENT"_together.tmp

bedtools sort -i "$SPECIES"_"$ELEMENT"_together.tmp > "$SPECIES"_"$ELEMENT"_together.sorted.tmp
 
bedtools merge -s -i "$SPECIES"_"$ELEMENT"_together.sorted.tmp -c 4 -o collapse > "$SPECIES"_"$ELEMENT"_together.merged.tmp

awk '{print $1 "\t" $2 "\t" $3 "\t" "BovB" "\t" "1" "\t" $4}' "$SPECIES"_"$ELEMENT"_together.merged.tmp > $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.bed

bedtools getfasta -s -fi $GENOME -bed $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.bed -fo $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.fasta

rm *tmp
