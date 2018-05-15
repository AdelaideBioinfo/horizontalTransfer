#!/bin/bash

# Example usage:
# GENOMEDIR=/data/rc003/atma/horizontalTransfer/genomes/fungi/Yarrowia.lipolytica GENOME=Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna QUERYDIR=/data/rc003/atma/horizontalTransfer/results/L1/tblastned QUERY=Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna_L1_ORFp.fasta_nucl_seqs.fasta RESULTSDIR=/data/rc003/atma/horizontalTransfer/results/L1/lastzed sbatch lastzAndExtract.q

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --time=0:10:00 
#SBATCH --mem=1GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au      

# Load BEDTools
module load BEDTools/2.25.0-foss-2015b

echo $GENOME
echo $QUERY

# Run lastz
lastz $GENOMEDIR/$GENOME[unmask,multiple] $QUERYDIR/$QUERY[unmask,multiple] --chain --gapped --coverage=80 --ambiguous=n --ambiguous=iupac --format=general-:name2,start2,end2,score,strand2,size2,name1,start1,end1 > $GENOME_$QUERY.lastzout

# Rearrange columns to put the concatenated file in BED-like form (for BEDTools) 
awk '{print $7 "\t" $8 "\t" $9 "\t" $1 "\t" "1" "\t" $5}' $GENOME_$QUERY.lastzout > BedFormat_$GENOME_$QUERY.tmp

# Sort by chr/scaffold and then by start position in ascending order 
bedtools sort -i BedFormat_$GENOME_$QUERY.tmp > Sorted_$GENOME_$QUERY.tmp

# Merge nested or overlapping intervals
bedtools merge -s -i Sorted_$GENOME_$QUERY.tmp -c 4 -o collapse > Merged_$GENOME_$QUERY.tmp

# Merging broke the BED-like format
# (i.e. strand is now in 4th column, merged names in 5th column, etc)
# Rearrange the columns as required
awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" "1" "\t" $4}' Merged_$GENOME_$QUERY.tmp > "$GENOME"_"$QUERY"_lastz.bed

# Extract FASTA from corrected merged BED file
bedtools getfasta -s -fi "$GENOMEDIR/$GENOME" -bed "$GENOME"_"$QUERY"_lastz.bed -fo "$GENOME"_"$QUERY"_lastz.fasta

# Move output bed and fasta to results folder
mv "$GENOME"_"$QUERY"_lastz.fasta $RESULTSDIR
mv "$GENOME"_"$QUERY"_lastz.bed $RESULTSDIR

# Remove temp files
rm *.tmp
rm *.lastzout
