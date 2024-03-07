
# Script to assign reads and barcodes to wt and mutant KRAS G12
# Daniele Tavernari (daniele.tavernari@epfl.ch)
# 06/03/2024

# Default KRAS G12 coordinates are for hg38
# Usage: inspect your samples' bam files on IGV, run the read extraction bash script first, adjust the 'Input' section and run this script.
# Output files:
# - *_KrasG12_MutationStatus_ReadLevel.txt: a table with read names, barcodes, reference codons, sequenced codons and mutation status for each read.
# - *_KrasG12_MutationStatus_BarcodeLevel.txt: a table with unique barcodes, reference codons, sequenced codons and mutation status for each barcode. Barcodes that are neither wt nor mutant are excluded. If a barcode had both a wt and a mutant read, the barcode is deemed mutant.
# The BarcodeLevel table can be read in downstream analysis to annotate cells into wt, mutant and unknown (i.e. cells with no reads overlapping KRAS G12, or neither wt nor mutant, and thus not appearing in the BarcodeLevel table).
# The ReadLevel table gives additional information into e.g. how much is the mutant form expressed.

######### Input
OutDir="/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/KRAS_G12_reads/"
InDir="/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/data/KRAS_G12_reads/"
SampleList="HRR338086 HRR338087 HRR338126 HRR338127 HRR338128 HRR338129 HRR338130 HRR338131 HRR338132 HRR338133 HRR338134 HRR338135 HRR338136 HRR338137"

g12_chr="chr12"
g12_start=25245349
g12_end=25245351

# See codon table here: https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables
wt_codons = c( "ACC","GCC","TCC","CCC" ) # as in IGV. Can be more than one (e.g. if you want to call synonymous variants 'wt'). 
mut_codons_1 = c( "ACA","GCA" ) # as in IGV. Can be more than one. 
mut_codons_2 = c( "ATC","GTC" )
mutant_name_1 = "G12C" # only two mutations are supported for now. If you want to accomodate for more mutations, you can modify this script or split mutations in downstream analyses.
mutant_name_2 = "G12D"
###############

SampleList = unlist(strsplit(SampleList,split=" "))
for (sample in SampleList){
	cat( "\n","Processing",sample,"...","\n" )
	fileName_barcodes = paste0( InDir,sample,"_barcodes.txt" )
	fileName_codons = paste0( InDir,sample,"_codons.txt" )
	if ((file.size(fileName_barcodes) > 0L) & (file.size(fileName_codons) > 0L) ){
		barcodes = read.table(fileName_barcodes,header=F,stringsAsFactors=F)
		codons = read.table(fileName_codons,header=F,stringsAsFactors=F,comment.char="")
		codons = codons[2:nrow(codons),]
		codons = codons[codons$V6!=".",c( "V1","V4","V8","V6","V9" )]
		colnames(codons) = c( "read_name","chr","pos","sequenced_base","reference_base" )
		colnames(barcodes) = c( "read_name","barcode" )
		barcodes$barcode = gsub( "CB:Z:","",barcodes$barcode )
		c1 = codons[codons$pos==g12_start,]
		c2 = codons[codons$pos==(g12_start+1),]
		c3 = codons[codons$pos==g12_end,]
		common_rn = intersect(intersect(barcodes$read_name,c1$read_name),intersect(c2$read_name,c3$read_name))
		barcodes = barcodes[barcodes$read_name %in% common_rn,]
		barcodes = barcodes[!duplicated(barcodes$read_name),]
		rownames(barcodes) = barcodes$read_name
		c1 = c1[c1$read_name %in% common_rn,]
		c1 = c1[!duplicated(c1$read_name),]
		rownames(c1) = c1$read_name
		c2 = c2[c2$read_name %in% common_rn,]
		c2 = c2[!duplicated(c2$read_name),]
		rownames(c2) = c2$read_name
		c3 = c3[c3$read_name %in% common_rn,]
		c3 = c3[!duplicated(c3$read_name),]
		rownames(c3) = c3$read_name
		barcodes[common_rn,"reference_codon"] = paste0(c1[common_rn,"reference_base"],c2[common_rn,"reference_base"],c3[common_rn,"reference_base"])
		barcodes[common_rn,"sequenced_codon"] = paste0(c1[common_rn,"sequenced_base"],c2[common_rn,"sequenced_base"],c3[common_rn,"sequenced_base"])
		barcodes$KRAS_G12 = "."
		barcodes[barcodes$sequenced_codon %in% wt_codons,"KRAS_G12"] = "wt"
		barcodes[barcodes$sequenced_codon %in% mut_codons_1,"KRAS_G12"] = mutant_name_1
		barcodes[barcodes$sequenced_codon %in% mut_codons_2,"KRAS_G12"] = mutant_name_2
		write.table( barcodes,file=paste0(OutDir,sample,"_KrasG12_MutationStatus_ReadLevel.txt"),sep="\t",quote=F,row.names=F,col.names=T )
		barcodes = barcodes[barcodes$KRAS_G12!=".",]
		barcodes$KRAS_G12 = factor(barcodes$KRAS_G12,levels = c( mutant_name_1,mutant_name_2,"wt" ))
		barcodes = barcodes[order(barcodes$KRAS_G12),]
		barcodes = barcodes[!duplicated(barcodes$barcode),]
		barcodes$read_name = NULL
		write.table( barcodes,file=paste0(OutDir,sample,"_KrasG12_MutationStatus_BarcodeLevel.txt"),sep="\t",quote=F,row.names=F,col.names=T )
	} else {
		cat( "\n","No reads with valid barcodes overlapping KRAS G12. No output file generated for",sample,"\n" )
	}
}




