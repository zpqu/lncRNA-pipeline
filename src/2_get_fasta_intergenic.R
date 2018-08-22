##This script is used to get fasta sequences for intergenic transcripts
##Note: Make sure there are reference genome in the DB folder
##Date: 17/02/2014
##Author: Zhipeng

library("Rsamtools")

###############function
myfun = function(file){
	fa = open(FaFile("../DB/genome.fa"))
	idx = scanFaIndex(fa)

	load(paste("../Rdata/", file, ".assembled.annotation.grl.RData",sep = ""))
	if(length(all.assembled.exon.gene.intergenic.grl) != 0){
	  all.assembled.exon.gene.intergenic.grl = endoapply(all.assembled.exon.gene.intergenic.grl, function(x) if(as.character(strand(x))[1] == "-"){rev(x)} else {x})
	}
	seqs = getSeq(fa, unlist(all.assembled.exon.gene.intergenic.grl, use.names = F), as.character = T)
	elt = rep(names(all.assembled.exon.gene.intergenic.grl), elementNROWS(all.assembled.exon.gene.intergenic.grl))
	seqs.fa = DNAStringSet(sapply(split(seqs, elt), paste, collapse = ""))
#	names(seqs.fa) = paste("colXcol", names(seqs.fa), sep = "_")
	writeXStringSet(seqs.fa, paste("../Output/", file, ".assembled.exon.gene.intergenic.fa", sep = ""))
}

setwd("../../Output")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}

