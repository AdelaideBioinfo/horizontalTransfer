#!/bin/bash

# Example usage:
# INDIR=/data/rc003/atma/horizontalTransfer/results/BovB/final FILE=Bos.taurus_BovB_final.fasta ID=80 sbatch usearchConsensus.q

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

# run usearch to generate consensus seqs
usearch -cluster_fast $FILE -threads 8 -qmask none -id 0."$ID" -uc "$FILE"_"$ID".uc -consout "$FILE"_"$ID".consensus

