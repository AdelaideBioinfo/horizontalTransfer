#!/bin/bash

# Example usage:
# INDIR=/data/rc003/atma/horizontalTransfer/results/BovB/final FILE=Bos.taurus_BovB_final.fasta MINSEQLENGTH=2000 MAXSEQLENGTH=4000 sbatch usearchSortByLength.q

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-00:00
#SBATCH --mem=32GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# go to input dir
cd $INDIR

# run usearch to sort by sequence length
# can also set min and max seq length cutoffs
usearch -sortbylength $FILE -minseqlength $MINSEQLENGTH -maxseqlength $MAXSEQLENGTH -fastaout $FILE.sorted

