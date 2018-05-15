#!/bin/bash

# Example usage:
#
# INDIR=/data/rc003/atma/horizontalTransfer/results/BovB/final FILE=allSpecies_BovB_centroids.fasta MINBLOCKSIZE=5 ALLOWEDGAPS=a sbatch alignRefineAndTree.q
#

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH --time=1-00:00 
#SBATCH --mem=32GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au 

# load module 
module load MUSCLE/3.8.31

# go to dir
cd $INDIR

# align sequences
# note: use -maxiters 2 if very large alignment
muscle -in $FILE -out ${FILE%.fasta}.afa

# Get blocks of conserved sequence from the alignment file
# Allowed gap positions is by default set to none
# to change, set -b5=h (half) or a (all)
# Min block size is by default set to 10

Gblocks ${FILE%.fasta}.afa -t=d -p=n -e=.gb -b4="$MINBLOCKSIZE" -b5="$ALLOWEDGAPS"
	
# Gblocks outputs an alignment file of 10-character blocks separated by white space
# Need to delete the white spaces
# Also change extension from .afa.gb to _gb.afa 
tr -d " \t" < ${FILE%.fasta}.afa.gb > ${FILE%.fasta}_b4"$MINBLOCKSIZE"_b5"$ALLOWEDGAPS"_gb.afa 

# Remove all files with white space and extension .afa.gb
#rm *.afa.gb

# Infer a maximum likelihood phylogeny from the refined alignment using FastTree
FastTree -nt -gtr < ${FILE%.fasta}_b4"$MINBLOCKSIZE"_b5"$ALLOWEDGAPS"_gb.afa > ${FILE%.fasta}_b4"$MINBLOCKSIZE"_b5"$ALLOWEDGAPS"_gb.tree 

