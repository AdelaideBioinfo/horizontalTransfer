#!/bin/bash

# Example usage:
# DIR=/data/rc003/atma/horizontalTransfer/genomes/fungi/Yarrowia.lipolytica DATABASE=Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna QUERY=L1_ORFp.fasta RESULTSDIR=/data/rc003/atma/horizontalTransfer/results/L1/tblastned sbatch tblastnAndExtract.q

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 2
#SBATCH --time=0:10:00 
#SBATCH --mem=1GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au      

# Load the necessary modules
module load BLAST+/2.6.0-foss-2016uofa-Python-2.7.11

tblastn -db "$DIR/$DATABASE" -num_threads 2 -query /data/rc003/atma/horizontalTransfer/queries/"$QUERY" -out "$DATABASE"_"$QUERY".out -outfmt "6 qseqid sseqid bitscore qcovs evalue pident length mismatch gapopen qstart qend qseq sstart send sseq" -evalue 1e-5 -show_gis \
2> "$DATABASE"_"$QUERY".log
	
cat "$DATABASE"_"$QUERY".out \
| sort -nr -k3 \
> "$DATABASE"_"$QUERY".out.sorted

# Put the target seq hits in FASTA format (for further checking)
# first, remove gaps ("-") from the target sequence (column 15)
# then print out all of the target seqs, with the header in the form: 
# >sseqid:sstart-send	qseqid:qstart-qend	score:*,qcov:*,evalue:*,pident:*
# e.g. >gi|699251403|ref|XM_002119377.3|:618-821	L1-1a_Cis_ORFp:208-291	score:38.5,qcov:17,evalue:6e-56,pid:30.68

cat "$DATABASE"_"$QUERY".out.sorted | awk -F" " '{gsub(/[-]/,"",$15)}1' OFS=" " | awk '{print ">" $2 ":" $13 "-" $14 "\t" $1 ":" $10 "-" $11 "\t" "score:" $3 ",qcov:" $4 ",eval:" $5 ",pid:" $6 "\n" $15}' | sed 's/ref|//g' | sed 's/|:/:/g' \
> "$DATABASE"_"$QUERY"_sseq.fasta

# Note: This output is the protein hit sequences, not the nucleotide sequences
# So now, extract nucleotide FASTA sequence for each hit

# Grab the headers from the amino acid FASTA file
grep '>' "$DATABASE"_"$QUERY"_sseq.fasta > "$DATABASE"_"$QUERY"_sseq_headers.txt

# Pull out all the plus strand sequences
cat "$DATABASE"_"$QUERY"_sseq_headers.txt \
| cut -f 1 \
| sed 's/>//g'  \
| sed 's/:/\t/g' \
| awk 'BEGIN {OFS=FS="\t"} {gsub(/\-/,"\t",$2)}1' \
| awk '{if ($2 < $3) print $1 "\t" $2 "\t" $3 "\t" "L1hit" "\t" "1" "\t" "+"}' \
> "$DATABASE"_"$QUERY"_plus.txt 

# Pull out all the minus strand sequences
cat "$DATABASE"_"$QUERY"_sseq_headers.txt \
| cut -f 1 \
| sed 's/>//g'  \
| sed 's/:/\t/g' \
| awk 'BEGIN {OFS=FS="\t"} {gsub(/\-/,"\t",$2)}1' \
| awk '{if ($2 > $3) print $1 "\t" $3 "\t" $2 "\t" "L1hit" "\t" "1" "\t" "-"}' \
> "$DATABASE"_"$QUERY"_minus.txt 

# Combine minus and plus strand sequences
cat "$DATABASE"_"$QUERY"_plus.txt "$DATABASE"_"$QUERY"_minus.txt  > "$DATABASE"_"$QUERY"_all.txt

# Load BEDTools
module load BEDTools/2.25.0-foss-2015b

# Sort by chr/scaffold and then by start position in ascending order 
bedtools sort -i "$DATABASE"_"$QUERY"_all.txt > "$DATABASE"_"$QUERY"_all.txt.sorted

# Merge nested or overlapping intervals
bedtools merge -s -i "$DATABASE"_"$QUERY"_all.txt.sorted -c 4 -o collapse > "$DATABASE"_"$QUERY"_all.txt.sorted.merged

# Merging broke the BED-like format
# (i.e. strand is now in 4th column, merged names in 5th column, etc)
# Rearrange the columns as required
awk '{print $1 "\t" $2 "\t" $3 "\t" "L1hit" "\t" "1" "\t" $4}' "$DATABASE"_"$QUERY"_all.txt.sorted.merged > "$DATABASE"_"$QUERY"_merged.bed

# Extract FASTA from corrected merged BED file
bedtools getfasta -s -fi "$DIR/$DATABASE" -bed "$DATABASE"_"$QUERY"_merged.bed -fo "$DATABASE"_"$QUERY"_results.fasta

# Sort sequences by length
usearch -sortbylength "$DATABASE"_"$QUERY"_results.fasta -fastaout "$DATABASE"_"$QUERY"_nucl_seqs.fasta

# Move output bed and fasta to results folder
mv "$DATABASE"_"$QUERY"_nucl_seqs.fasta $RESULTSDIR
mv "$DATABASE"_"$QUERY"_merged.bed $RESULTSDIR

# Remove or move remaining files to a temporary folder
rm "$DATABASE"_"$QUERY"* 
