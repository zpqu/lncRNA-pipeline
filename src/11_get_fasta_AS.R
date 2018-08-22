##This script is used to get fasta sequences for intergenic transcripts
##Date: 12/05/2014
##Author: Zhipeng

library("Rsamtools")

###############function
myfun = function(file){
	fa = open(FaFile("../DB/genome.fa"))
	idx = scanFaIndex(fa)

	load(paste("../Rdata/", file, ".assembled.annotation.grl.RData",sep = ""))
	if(length(all.assembled.exon.gene.AS.grl) != 0){
        all.assembled.exon.gene.AS.grl = endoapply(all.assembled.exon.gene.AS.grl, function(x) if(as.character(strand(x))[1] == "-"){rev(x)} else {x})
	}
	if(length(all.assembled.exon.gene.intronic.AS.grl) != 0){
	all.assembled.exon.gene.intronic.AS.grl = endoapply(all.assembled.exon.gene.intronic.AS.grl, function(x) if(as.character(strand(x))[1] == "-"){rev(x)} else {x})
	}
	if(length(all.assembled.exon.gene.intronic.grl) != 0){
	all.assembled.exon.gene.intronic.grl = endoapply(all.assembled.exon.gene.intronic.grl, function(x) if(as.character(strand(x))[1] == "-"){rev(x)} else {x})
	}

	####For AS_exon
	seqs = getSeq(fa, unlist(all.assembled.exon.gene.AS.grl, use.names = F), as.character = T)
	elt = rep(names(all.assembled.exon.gene.AS.grl), elementNROWS(all.assembled.exon.gene.AS.grl))
	seqs.fa = DNAStringSet(sapply(split(seqs, elt), paste, collapse = ""))
#	names(seqs.fa) = paste("colXcol", names(seqs.fa), sep = "_")
	writeXStringSet(seqs.fa, paste("../Output/", file, ".assembled.exon.gene.AS.fa", sep = ""))

	####For AS_intronic
	seqs = getSeq(fa, unlist(all.assembled.exon.gene.intronic.AS.grl, use.names = F), as.character = T)
	elt = rep(names(all.assembled.exon.gene.intronic.AS.grl), elementNROWS(all.assembled.exon.gene.intronic.AS.grl))
	seqs.fa = DNAStringSet(sapply(split(seqs, elt), paste, collapse = ""))
#	names(seqs.fa) = paste("colXcol", names(seqs.fa), sep = "_")
	writeXStringSet(seqs.fa, paste("../Output/", file, ".assembled.exon.gene.intronic.AS.fa", sep = ""))

	###For intronic
	seqs = getSeq(fa, unlist(all.assembled.exon.gene.intronic.grl, use.names = F), as.character = T)
	elt = rep(names(all.assembled.exon.gene.intronic.grl), elementNROWS(all.assembled.exon.gene.intronic.grl))
	seqs.fa = DNAStringSet(sapply(split(seqs, elt), paste, collapse = ""))
#	names(seqs.fa) = paste("colXcol", names(seqs.fa), sep = "_")
	writeXStringSet(seqs.fa, paste("../Output/", file, ".assembled.exon.gene.intronic.fa", sep = ""))


}

setwd("../../Output")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf", "", files[i])
      myfun(file = fname)
}
