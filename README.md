# Eukaryotic Horizontal Transfer

Code used to infer and filter horizontal transfer events involving L1 and BovB retrotransposons in eukaryotes. 

Also available on Zenodo:

- Atma M Ivancevic. (2018, May 15). AdelaideBioinfo/horizontalTransfer: First release of horizontal transfer code (Version v1.0.0). Zenodo. http://doi.org/10.5281/zenodo.1246999

- Ivancevic, Atma M, Kortschak, R Daniel, Bertozzi, Terry, & Adelson, David L. (2018). Dataset from: Horizontal transfer of BovB and L1 retrotransposons in eukaryotes [Data set]. Zenodo. http://doi.org/10.5281/zenodo.1246946

#### Recommended programs
- BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- SAMtools (http://www.htslib.org/)
- BEDtools (http://bedtools.readthedocs.io/en/latest/)
- CENSOR, which requires wu-blast and bioperl (https://girinst.org/downloads/software/censor/)
- USEARCH (https://www.drive5.com/usearch/)
- VSEARCH (https://github.com/torognes/vsearch)
- MUSCLE (https://www.drive5.com/muscle/)

#### Optional
- LASTZ (https://www.bx.psu.edu/~rsharris/lastz/)
- SiLiX (http://lbbe.univ-lyon1.fr/-SiLiX-?lang=en)
- RepeatMasker (http://www.repeatmasker.org/)
- Gblocks (https://ls23l.lscore.ucla.edu/MakeTree/documentation/gblocks.html)
- HMMer (http://hmmer.org/)
- FastTree (http://www.microbesonline.org/fasttree/)

#### Prerequisites
- Some level of familiarity with queuing systems (SLURM)

#### Scripts and test genome
A test genome (fungus *Yarrowia lipolytica*) has been placed in [genomes/fungi/Yarrowia.lipolytica](genomes/fungi/Yarrowia.lipolytica). We recommend trying out the workflow below on this genome first. Intermediate files for each step can be found in [results/L1](results/L1), to help with troubleshooting. Analysis scripts are provided in [scripts](scripts). LaTeX files for the Additional Files (Supp Tables and Figures) are provided in [latex](latex).

## Workflow

### 1) Extraction of L1 and BovB retrotransposons from genome data

#### 1a) Download or acquire genomes of interest
755 genomes were downloaded from public databases (UCSC and NCBI); four more were acquired from collaborators. All genomes were downloaded using ```wget```, as recommended by NCBI (https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/). See Supplementary Table 1 for the source and assembly version of each genome used.

#### 1b) Append species name to each genome fasta file
This will save you a lot of headaches later on. Most genomes are given useless names like GCF_000002525.2_ASM252v1_genomic.fna.
Use [renameGenome.q](scripts/renameGenome.q). 

For example, to rename fungal genome *Yarrowia lipolytica*:
```
GROUP=fungi SPECIES=Yarrowia.lipolytica GENOME=GCF_000002525.2_ASM252v1_genomic.fna sbatch renameGenome.q
```

#### 1c) Make each genome a BLAST database and create index file
See [makeDatabaseAndIndex.q](scripts/makeDatabaseAndIndex.q). General commands are:
```
makeblastdb -in genome.fna -parse_seqids -dbtype nucl

samtools faidx genome.fna
```

#### 1d) Find TEs using TBLASTN with protein queries
[tblastnAndExtract.q](scripts/tblastnAndExtract.q) will run tblastn to find protein matches, then extract the corresponding nucleotide sequences from the genome. 
- DIR is the directory where your genome is stored.
- DATABASE is the name of your genome file (which we used to create a BLAST database).
- QUERY is the file containing the protein sequences you're searching for (e.g. [L1_ORFp.fasta](queries/L1_ORFp.fasta) or [BovB_ORFp.fasta](queries/BovB_ORFp.fasta)).
- RESULTSDIR is your specified output directory.

For example, to run [tblastnAndExtract.q](scripts/tblastnAndExtract.q) on fungal genome *Yarrowia lipolytica* with L1 queries:
```
DIR=/data/rc003/atma/horizontalTransfer/genomes/fungi/Yarrowia.lipolytica DATABASE=Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna QUERY=L1_ORFp.fasta RESULTSDIR=/data/rc003/atma/horizontalTransfer/results/L1/tblastned sbatch tblastnAndExtract.q
```
Adjust memory and cores as necessary for each genome. Mammalian genomes need a bit more oomph, use:
```
#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH --time=3-00:00 
#SBATCH --mem=32GB 
```

#### 1e) OPTIONAL Find TEs using LASTZ with nucleotide queries
[lastzAndExtract.q](scripts/lastzAndExtract.q) will run LASTZ to find nucleotide matches, then extract the corresponding nucleotide TE (L1 or BovB) sequences from the genome. 
- GENOMEDIR is the directory where your genome is stored.
- GENOME is the name of your genome
- QUERYDIR is the directory containing your query file (e.g. L1 nucleotide sequences)
- QUERY is your query file (e.g. L1 nucleotide sequences)
- RESULTSDIR is your specified output directory

E.g. to run [lastzAndExtract.q](scripts/lastzAndExtract.q) on fungal genome *Yarrowia lipolytica* using the L1 queries found above with TBLASTN:
```
GENOMEDIR=/data/rc003/atma/horizontalTransfer/genomes/fungi/Yarrowia.lipolytica GENOME=Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna QUERYDIR=/data/rc003/atma/horizontalTransfer/results/L1/tblastned QUERY=Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna_L1_ORFp.fasta_nucl_seqs.fasta RESULTSDIR=/data/rc003/atma/horizontalTransfer/results/L1/lastzed sbatch lastzAndExtract.q
```
Large genomes (e.g. mammals) need to be split into smaller files before being passed through LASTZ.

**Note: we will probably stop using LASTZ for future analyses because we have not noticed a significant increase in true TE hits, even when used iteratively. This step is optional.** 

#### 1f) Combine LASTZ and TBLASTN hits
Use [combineHits.q](scripts/combineHits.q) to combine LASTZ and TBLASTN hits, for further checking. 

E.g. to combine all L1 hits for genome *Yarrowia lipolytica*:
```
SPECIES=Yarrowia.lipolytica ELEMENT=L1 LASTZFILE=/data/rc003/atma/horizontalTransfer/results/L1/lastzed/Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna_Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna_L1_ORFp.fasta_nucl_seqs.fasta_lastz.bed TBLASTNFILE=/data/rc003/atma/horizontalTransfer/results/L1/tblastned/Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna_L1_ORFp.fasta_merged.bed GENOME=/data/rc003/atma/horizontalTransfer/genomes/fungi/Yarrowia.lipolytica/Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna RESULTSDIR=/data/rc003/atma/horizontalTransfer/results/L1/combined sbatch combineHits.q
```

#### 1g) Use CENSOR to categorise all extracted nucleotide sequences
[censorSequences.q](scripts/censorSequences.q) will use CENSOR to check sequences against the RepBase library of known repeat sequences. We want to see if our "L1 hits" are truly L1s, as deemed by CENSOR, and if our "BovB hits" are truly BovBs. Sometimes, other repeats are picked up instead.
- DIR is the directory where you put your combined tblastn+lastz results
- FILE is the combined fasta file
- RESULTSDIR is the directory that you want the CENSOR output to go in

E.g. to run the CENSOR on all L1 hits from fungal genome *Yarrowia lipolytica*:
```
DIR=/data/rc003/atma/horizontalTransfer/results/L1/combined FILE=Yarrowia.lipolytica_L1_combined.fasta RESULTSDIR=/data/rc003/atma/horizontalTransfer/results/L1/censored sbatch censorSequences.q
```

#### 1h) Confirm hits are really L1 or BovB elements
[checkCensorOutput.q](scripts/checkCensorOutput.q) uses the output .map file from CENSOR to find which hits are recognised as L1 or BovB sequences. 
- SPECIES is the species name
- FILE is the tblastn results file
- GENOME is the genome you are querying
- ELEMENT is your retrotransposon (ie. L1 or BovB)
- QUERY is the text file of RepBase elements you are using to confirm your L1/BovB hits (e.g. [1826_L1_and_Tx1_repbase_headers.txt](queries/1826_L1_and_Tx1_repbase_headers.txt) or [37_BovB_headers.txt](queries/37_BovB_headers.txt))

E.g. to confirm L1 hits in fungal genome *Yarrowia lipolytica*:
```
SPECIES=Yarrowia.lipolytica FILE=Yarrowia.lipolytica_L1_combined.fasta GENOME=/data/rc003/atma/horizontalTransfer/genomes/fungi/Yarrowia.lipolytica/Yarrowia.lipolytica_GCF_000002525.2_ASM252v1_genomic.fna ELEMENT=L1 QUERY=/data/rc003/atma/horizontalTransfer/queries/1826_L1_and_Tx1_repbase_headers.txt sbatch checkCensorOutput.q
```

#### 1i) Make all fasta uppercase
Some programs automatically mask or ignore lowercase nucleotides (i.e. repetitive regions). To avoid bizarre clusterings later, convert all L1 and BovB sequences to uppercase using [toUpper.q](scripts/toUpper.q).

E.g. to upperize L1 hits in fungal genome *Yarrowia lipolytica*:
```
ELEMENT=L1 SPECIES=Yarrowia.lipolytica sbatch toUpper.q
```

#### 1j) Append species name to sequence headers
So far we have dealing with each genome separately. Later, we will combine ALL detected L1 sequences from all genomes into one monster file. This will allow us to perform all-against-all clustering to find discordant groupings.

Appending the species name to each sequence header makes it easier to keep track of which species/genome each TE came from. Use [appendNameToHeaders.q](scripts/appendNameToHeaders.q), 

e.g. to append the species name "Yarrowia.lipolytica" to all L1 sequences from *Yarrowia lipolytica*:
```
ELEMENT=L1 SPECIES=Yarrowia.lipolytica sbatch appendNameToHeaders.q
```

This script also moves all final uppercased, renamed L1 sequences to the directory `final/`. 

This is the end of Part 1. Repeat from the beginning to find and extract BovB sequences (i.e. change the QUERY and ELEMENT variables to BovB).

### 2) Inferring a representative BovB tree from consensus/centroid sequences

#### 2a) Consensus sequence approach
[usearchConsensus.q](scripts/usearchConsensus.q) was used to generate consensus sequence(s) for each species containing BovB. 

E.g. to generate consensus sequences from the final BovB seqs in cow *Bos taurus*:
```
INDIR=/data/rc003/atma/horizontalTransfer/results/BovB/final FILE=Bos.taurus_BovB_final.fasta ID=80 sbatch usearchConsensus.q
```

#### 2b) Finding the longest BovB per species
[usearchSortByLength.q](scripts/usearchSortByLength.q) was used to sort the BovBs in each genome from largest to smallest. If the longest BovB was ridiculously long (e.g. 20kb, from many duplicated copies next to each other), then a maximum seq length of 4kb was adopted. Min seq length for BovBs was set to 2kb.

E.g. *Bos taurus* BovBs were sorted by length, with min and max length cutoffs of 2kb and 4kb respectively:
```
INDIR=/data/rc003/atma/horizontalTransfer/results/BovB/final FILE=Bos.taurus_BovB_final.fasta MINSEQLENGTH=2000 MAXSEQLENGTH=4000 sbatch usearchSortByLength.q
```
The longest BovB under 4kb was chosen as the species "centroid". Centroids were then concatenated into one file:
```
cat *_BovB_final.fasta.centroid > allSpecies_BovB_centroids.fasta
```

#### 2c) Align, refine and infer tree
[alignRefineAndTree.q](scripts/alignRefineAndTree.q) was used to generate a MUSCLE alignment of the BovB centroids, use Gblocks to extract conserved alignment blocks, and infer a maximum likelihood tree with FastTree:
```
INDIR=/data/rc003/atma/horizontalTransfer/results/BovB/final FILE=allSpecies_BovB_centroids.fasta MINBLOCKSIZE=5 ALLOWEDGAPS=a sbatch alignRefineAndTree.q
```

### 3) Distinguishing between RTE and BovB elements

For each genome, used [usearchCluster.q](scripts/usearchCluster.q) to cluster all BovB and other identified RTE sequences together. This determined which "RTE" sequences from RepBase were actually BovBs.

### 4) Clustering of nucleotide BovB sequences from bats and frogs

As above, concatenated all frog and bat BovB sequences; then used [usearchSortByLength.q](scripts/usearchSortByLength.q) to sort by sequence length (min 2kb, max 4kb); then used [usearchCluster.q](scripts/usearchCluster.q) to cluster the sequences; and [alignRefineAndTree.q](scripts/alignRefineAndTree.q) to align and infer a phylogeny. 

### 5) Extraction of nucleotide ORFs and conserved amino acid residues

#### 5a) Extract open reading frames and reverse-transcriptase domains

[findOrfs.q](scripts/findOrfs.q) was used to find open reading frames, and extract:
- nucleotide open reading frames
- amino acid open reading frames
- amino acid reverse transcriptase domains

ORFs and RT domains can be extracted per species, or from the concatenated L1/BovB multi-genome files. 
E.g. to find ORFs present in our L1 dataset, with min codons 600:
```
INDIR=/data/rc003/atma/horizontalTransfer/results/L1/allSpeciesTogether FILE=allSpecies_L1.fasta MINCODONS=600 sbatch findOrfs.sh
```
#### 5b) OPTIONAL clustering of ORFs and RT domains
These ORF and RT files were then clustered, to find potential instances of horizontal transfer.
[vsearchCluster.q](scripts/vsearchCluster.q) was used for ORF nucleotide sequences.
[usearchCluster.q](scripts/usearchCluster.q) was used for RT amino acid sequences.
See the clustering section below for more details.

**Important note: vsearch does not currently support protein sequences. However, it will not fail if you give it a protein sequence input - you will just get a really terrible clustering as it looks for nucleotide bases.**

### 6) All-against-all clustering

#### 6a) Make separate multi-fasta databases for L1 and BovB elements
L1s:
```
# go to dir
cd /data/rc003/atma/horizontalTransfer/results/L1/final
# concatenate L1s from all genomes
cat *_L1_final.fasta > ../allSpeciesTogether/allSpecies_L1.fasta
# sort and output sequences between 3-9kb
vsearch -sortbylength allSpecies_L1.fasta -minseqlength 3000 -maxseqlength 9000 --output allSpecies_L1_min3kb_max9kb.fasta
```

Likewise for BovBs:
```
# go to dir
cd /data/rc003/atma/horizontalTransfer/results/BovB/final
# concatenate BovBs from all genomes
cat *_BovB_final.fasta > ../allSpeciesTogether/allSpecies_BovB.fasta
# sort and output sequences between 2.4-4kb
vsearch -sortbylength allSpecies_BovB.fasta -minseqlength 2400 -maxseqlength 4000 --output allSpecies_BovB_min2.4kb_max4kb.fasta
```

#### 6b) Cluster away!
Use [vsearchCluster.q](scripts/vsearchCluster.q) to cluster each database at different identities. 

E.g. to perform an all-against-all clustering of the L1 database at identity 80%, with cluster prefix "c_":
```
INDIR=/data/rc003/atma/horizontalTransfer/results/L1/allSpeciesTogether FILE=allSpecies_L1_min3kb_max9kb.fasta ID=80 PREFIX=c sbatch vsearchCluster.q
```
For the paper, we trialled clustering identities of 50-90%. We used "c" as the prefix for nucleotide sequence clusters, "o" as the prefix for open reading frame clusters, and "r" as the prefix for reverse-transcriptase domain clusters. 

We also tried using BLAST+SiLiX for the all-against-all clustering, as described previously (https://github.com/AdelaideBioinfo/L1-dynamics). However, this proved to be orders of magnitude slower with no apparent improvement over VSEARCH.

#### 6c) Filter clusters

Use [filterClusters.q](scripts/filterClusters.q) to categorize clusters based on whether they represent HT events which are cross-Species, cross-Order, cross-Class or cross-Phylum. The cross-Class and cross-Phylum ones tend to be the most interesting.

E.g. run the following after [vsearchCluster.q](scripts/vsearchCluster.q) to categorise L1 nucleotide seq clusters:
```
INDIR=/data/rc003/atma/horizontalTransfer/results/L1/allSpeciesTogether FILE=allSpecies_L1_min3kb_max9kb.fasta ID=80 PREFIX=c sbatch filterClusters.q
```

### 7) Kimura divergence estimates for species containing both TEs

BovB and L1 sequences for each genome were concatenated and screened with RepeatMasker using [repeatmasker.q](scripts/repeatmasker.q).

E.g. To run RepeatMasker on cow *Bos taurus*:
```
FILE=Bos.taurus_bothL1andBovB.fasta sbatch repeatmasker.q
```

Use [plotKimuraDivergence.R](scripts/plotKimuraDivergence.R) to plot the corresponding Kimura divergence of L1s and BovBs.

## THE END

In the words of George Box,
> "Essentially, all models are wrong, but some are useful"
