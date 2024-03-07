# KrasG12_MutCaller_10xReads

Scripts to identify and extract KRAS G12 mutations from 10x scRNA-seq reads. 

## 1. Preliminary steps
* inspect your samples' bam files on IGV to make sure you see the mutation, and that the coordinates are correct
* download [jvarkit](https://github.com/lindenb/jvarkit) (precompiled .jar available under https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH? ). Credits to Pierre Lindenbaum, PhD.
* make sure you have a recent java version. In case, create a conda environment, activate it (conda activate jvarkit) and install a new one
* create a sequence dictionary for your cellranger reference genome (you will need a recent samtools version) using samtools dict

## 2. Run extract_KrasG12_reads.sh
See usage instructions in the script

## 3. Run assignBarcodes_KrasG12mutant.R
See usage instructions in the script
