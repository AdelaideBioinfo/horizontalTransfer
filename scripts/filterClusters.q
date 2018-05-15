#!/bin/bash

# Example usage:
# INDIR=/data/rc003/atma/horizontalTransfer/results/L1/allSpeciesTogether FILE=allSpecies_L1_min3kb_max9kb.fasta ID=80 PREFIX=c sbatch filterClusters.q

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8 
#SBATCH --time=3-00:00 
#SBATCH --mem=30GB 

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

echo $(date +"[%b %d %H:%M:%S] Start")

cd $INDIR/
cd "$FILE"_"$ID"_clusters/

# summarise each cluster
for i in "$PREFIX"_*; do cat $i | grep '>' | sed 's/>//g' | awk -F "_" '{print $1}' | sort | uniq | awk '{print "'$i'" " " $1}' > species_"$i"; done

echo $(date +"[%b %d %H:%M:%S] Summarised clusters")

# find which cluster summary files contain >1 line (indicating more than one species)
wc -l species_"$PREFIX"_* | awk '{if ($1!=1) print}' | awk '{print $2}' > moreThanOne.txt

echo $(date +"[%b %d %H:%M:%S] Made moreThanOne.txt")

# keep these and move to a new dir
mkdir -p ../id"$ID"_moreThanOneSpecies_"$FILE"_clusters

while IFS='' read -r line || [[ -n "$line" ]]; do mv $line ../id"$ID"_moreThanOneSpecies_"$FILE"_clusters/; done < moreThanOne.txt

cd ../id"$ID"_moreThanOneSpecies_"$FILE"_clusters/

echo $(date +"[%b %d %H:%M:%S] Moved to Species dir")

##################################
#Species to Order

# replace species names with order (since we're restricting possible HT events to those that cross-Order
for i in species_"$PREFIX"_*;
do
awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1' /data/rc003/atma/horizontalTransfer/references/lists/speciesToOrderDict.txt "$i" > "$i"_order.tmp;
done

echo $(date +"[%b %d %H:%M:%S] Replaced Species with Order")

# remove duplicates to condense it down just to the orders present
for i in species_"$PREFIX"_*_order.tmp; do cat $i | sort | uniq > "${i%.tmp}"; done

echo $(date +"[%b %d %H:%M:%S] Removed duplicates")

# find clusters which have at least two different Orders
wc -l species_"$PREFIX"_*_order | awk '{if ($1!=1) print}' | awk '{print $2}' > moreThanOneOrder.txt

echo $(date +"[%b %d %H:%M:%S] Made moreThanOneOrder.txt")

# keep these and move to a new dir
mkdir -p ../id"$ID"_moreThanOneOrder_"$FILE"_clusters/

while IFS='' read -r line || [[ -n "$line" ]]; do mv $line ../id"$ID"_moreThanOneOrder_"$FILE"_clusters/; done < moreThanOneOrder.txt

cd ../id"$ID"_moreThanOneOrder_"$FILE"_clusters/

echo $(date +"[%b %d %H:%M:%S] Moved to Order dir")

##################################
#Order to Class

# replace names 
for i in species_"$PREFIX"_*_order;
do
awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1' /data/rc003/atma/horizontalTransfer/references/lists/orderToClassDict.txt "$i" > "$i"_class.tmp;
done

echo $(date +"[%b %d %H:%M:%S] Replaced Order with Class")

# remove duplicates 
for i in species_"$PREFIX"_*_class.tmp; do cat $i | sort | uniq > "${i%.tmp}"; done

echo $(date +"[%b %d %H:%M:%S] Removed duplicates")

# find clusters which have at least two different combinations
wc -l species_"$PREFIX"_*_class | awk '{if ($1!=1) print}' | awk '{print $2}' > moreThanOneClass.txt

echo $(date +"[%b %d %H:%M:%S] Made moreThanOneClass.txt")

# keep these and move to a new dir
mkdir -p ../id"$ID"_moreThanOneClass_"$FILE"_clusters

while IFS='' read -r line || [[ -n "$line" ]]; do mv $line ../id"$ID"_moreThanOneClass_"$FILE"_clusters/; done < moreThanOneClass.txt

cd ../id"$ID"_moreThanOneClass_"$FILE"_clusters/

echo $(date +"[%b %d %H:%M:%S] Moved to Class dir")

##################################
#Class to Phylum

# replace names 
for i in species_"$PREFIX"_*_class;
do
awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1' /data/rc003/atma/horizontalTransfer/references/lists/classToPhylumDict.txt "$i" > "$i"_phylum.tmp;
done

echo $(date +"[%b %d %H:%M:%S] Replaced Class with Phylum")

# remove duplicates 
for i in species_"$PREFIX"_*_phylum.tmp; do cat $i | sort | uniq > "${i%.tmp}"; done

echo $(date +"[%b %d %H:%M:%S] Removed duplicates")

# find clusters which have at least two different combinations
wc -l species_"$PREFIX"_*_phylum | awk '{if ($1!=1) print}' | awk '{print $2}' > moreThanOnePhylum.txt

echo $(date +"[%b %d %H:%M:%S] Made moreThanOnePhylum.txt")

# keep these and move to a new dir
mkdir -p ../id"$ID"_moreThanOnePhylum_"$FILE"_clusters

while IFS='' read -r line || [[ -n "$line" ]]; do mv $line ../id"$ID"_moreThanOnePhylum_"$FILE"_clusters/; done < moreThanOnePhylum.txt

cd ../id"$ID"_moreThanOnePhylum_"$FILE"_clusters/

echo $(date +"[%b %d %H:%M:%S] This is the end")
