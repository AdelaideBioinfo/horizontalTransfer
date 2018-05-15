#!/bin/bash

# Example usage:
# INDIR=/data/rc003/atma/horizontalTransfer/results/L1/allSpeciesTogether FILE=allSpecies_L1_min3kb_max9kb.fasta ID=50 PREFIX=c sbatch vsearchCluster.q

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1-00:00
#SBATCH --mem=32GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# go to input dir
cd $INDIR

# make dir for clusters
mkdir -p "$FILE"_"$ID"_clusters

# run vsearch
vsearch -cluster_fast $FILE -threads 16 -qmask none -id 0."$ID" -uc "$FILE"_"$ID".uc -clusters "$FILE"_"$ID"_clusters/"$PREFIX"_
