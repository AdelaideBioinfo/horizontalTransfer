#!/bin/bash

# Example usage:
# FILE=Bos.taurus_bothL1andBovB.fasta sbatch repeatmasker.q

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8 
#SBATCH --time=3-00:00 
#SBATCH --mem=32GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au 

# Run repeatmasker
RepeatMasker -pa 8 -a -qq -nolow -species eukaryotes $FILE

# Calculate the divergence from alignment output
# -s is the output of divergence test
perl /data/rc003/atma/RepeatMasker/util/calcDivergenceFromAlign.pl -s $FILE.divsum $FILE.align

