#!/bin/bash

# Rename genome file to include species name
#
# Example usage:
# GROUP=fungi SPECIES=Yarrowia.lipolytica GENOME=GCF_000002525.2_ASM252v1_genomic.fna sbatch renameGenome.q

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

cd /data/rc003/atma/horizontalTransfer/genomes/$GROUP/$SPECIES

mv $GENOME "$SPECIES"_"$GENOME"
