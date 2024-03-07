
# Script to extract 10x reads overlapping KRAS G12, their actual codons and their valid barcodes
# Daniele Tavernari (daniele.tavernari@epfl.ch)
# 06/03/2024

# Default KRAS G12 coordinates are for hg38
# Usage: follow the 'Preliminary steps', adjust the 'Input' section and run the script
# Output files:
# - *_codons.txt: a table with read names, sequenced and reference base for the 3 coordinates of the input codon. Note that many reads will not have a sequenced base for the codon coordinates, because they overlap it with their gap (see IGV).
# - *_barcodes.txt: a table with read names and their corresponding cell barcodes. Note that only a subset of reads (even zero) have a valid barcode (the others will be discarded anyway by cellranger when quantifying gene expression. See https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-bam)

######### Preliminary steps
# 1. inspect your samples' bam files on IGV to make sure you see the mutation, and that the coordinates are correct
# 2. download jvarkit precompiled .jar from https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH?
# 3. make sure you have a recent java version. In case, create a conda environment, activate it (conda activate jvarkit) and install a new one
# 4. create a sequence dictionary for your cellranger reference genome (you will need a recent samtools version): ${samtools_path}samtools dict ${refGenome_path}"genome.fa" -o ${refGenome_path}"genome.dict"
###########################

######### Input
OutDir="/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/KRAS_G12_reads/"
InDir="/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/Wang2021_HRA001130/cellranger_redo/"
SampleList="HRR338086 HRR338087 HRR338126 HRR338127 HRR338128 HRR338129 HRR338130 HRR338131 HRR338132 HRR338133 HRR338134 HRR338135 HRR338136 HRR338137"

g12_chr="chr12"
g12_start="25245349"
g12_end="25245351"

samtools_path="/mnt/ndata/daniele/lung_multiregion/Scripts/Tools/samtools-1.9/"
jvarkit_path="/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/utils/"
refGenome_path="/mnt/ndata/daniele/lung_multiregion/Scripts/Tools/refdata-gex-GRCh38-2020-A/fasta/"
###############

mkdir -p ${OutDir}
printf ${g12_chr}'\t'${g12_start}'\t'${g12_end}'\n' > ${OutDir}"g12_coordinates.bed"

for sample in ${SampleList}
do
	echo ${sample}
	${samtools_path}samtools view -b ${InDir}${sample}"/outs/possorted_genome_bam.bam" ${g12_chr}':'${g12_start}'-'${g12_end} > ${OutDir}"temp.bam"
	${samtools_path}samtools index ${OutDir}"temp.bam"
	${samtools_path}samtools view ${OutDir}"temp.bam" | cut -f 1,25 | grep "CB:Z" > ${OutDir}${sample}"_barcodes.txt"
	java -jar ${jvarkit_path}jvarkit.jar sam2tsv --reference ${refGenome_path}"genome.fa" --regions ${OutDir}"g12_coordinates.bed" -o ${OutDir}${sample}"_codons.txt" ${OutDir}"temp.bam"
	awk -v min=${g12_start} -v max=${g12_end} '(NR == 1) || ($8 >= min && $8 <= max) { print }' ${OutDir}${sample}"_codons.txt" > ${OutDir}${sample}"_codons.temp" && mv ${OutDir}${sample}"_codons.temp" ${OutDir}${sample}"_codons.txt"
	rm ${OutDir}temp.bam*
done

