#!/bin/bash

# Example usage:
# DIR=/data/rc003/atma/horizontalTransfer/results/L1/combined FILE=Yarrowia.lipolytica_L1_combined.fasta RESULTSDIR=/data/rc003/atma/horizontalTransfer/results/L1/censored sbatch censorSequences.q

#SBATCH -A robinson
#SBATCH --qos=long
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=7-00:00
#SBATCH --mem=32GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au      

# Load the necessary modules
module load bioperl/1.6.924
module load censor/4.2.29
module load wu-blast/2.0
module load BEDTools/2.25.0-foss-2015b

# Run CENSOR on the extracted nucleotide seqs, using the Repbase library
# Note: download the latest Repbase library from http://www.girinst.org/repbase/
censor -bprm cpus=16 -lib /data/rc003/atma/horizontalTransfer/references/repBaseLibrary/RepBase21.03_all_seqs.ref "$DIR/$FILE"

# Move .map file to results directory
mv "$FILE".map $RESULTSDIR

# Remove or move temporary files
rm "$FILE".*
rm censor.*.log
