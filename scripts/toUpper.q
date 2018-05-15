#!/bin/bash

# Convert all lower case to upper case (ignoring seq headers)
#
# Example usage:
# ELEMENT=L1 SPECIES=Yarrowia.lipolytica sbatch toUpper.q

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=0:10:00
#SBATCH --mem=1GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au 

cd /data/rc003/atma/horizontalTransfer/results/$ELEMENT/censored

cat "$SPECIES"_"$ELEMENT"_censored.fasta \
| awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' \
> "$SPECIES"_"$ELEMENT"_upper.fasta
