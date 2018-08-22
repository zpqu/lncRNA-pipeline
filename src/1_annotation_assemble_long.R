##This script is for annotating reconstructed transcripts.
##Date:14/02/2014
##Author:Zhipeng Qu

rm(list = ls())
library("rtracklayer")
library("GenomicRanges")

print(paste("The program starts at: ", Sys.time()))
lincRNA.path = "../.."
load(paste(lincRNA.path, "/DB/geneModel.RData", sep = ""))

myfun = function(fname){
   
print(paste("Now is processing: ", fname, Sys.time()))
filename = gsub("\\.gtf", "", fname)
all.assembled = import(paste(lincRNA.path, "/GTF/", fname, sep = ""))
#seqlevels(all.assembled) = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
#seqnames(all.assembled) = gsub("chloroplast", "ChrC", seqnames(all.assembled))
#seqnames(all.assembled) = gsub("mitochondria", "ChrM", seqnames(all.assembled))
all.chr = gsub("chloroplast", "ChrC", as.character(seqnames(all.assembled)))
all.chr = gsub("mitochondria", "ChrM", all.chr)
all.assembled = GRanges(seqnames = Rle(all.chr), ranges = ranges(all.assembled), strand = strand(all.assembled), values(all.assembled))

all.assembled.exon.gr = all.assembled[elementMetadata(all.assembled)$type == "exon"]
all.assembled.exon.grl = split(all.assembled.exon.gr, elementMetadata(all.assembled.exon.gr)$transcript_id)

all.assembled.exon.length = sum(width(all.assembled.exon.grl))
save(all.assembled.exon.length, file = paste(lincRNA.path, "/Rdata/", filename, ".assembled.annotation.exon.length.RData", sep = ""))
all.assembled.exon.all.grl = all.assembled.exon.grl
all.assembled.exon.grl = all.assembled.exon.grl[all.assembled.exon.length > 200]

all.assembled.exon.gene.ol = findOverlaps(all.assembled.exon.grl, geneModel.exon.grl)
length(unique(queryHits(all.assembled.exon.gene.ol)))
length(unique(subjectHits(all.assembled.exon.gene.ol)))
all.assembled.exon.gene.grl = all.assembled.exon.grl[unique(queryHits(all.assembled.exon.gene.ol))]
all.assembled.exon.geneRemove.grl = subset(all.assembled.exon.grl, !(names(all.assembled.exon.grl) %in% names(all.assembled.exon.gene.grl)))

all.assembled.exon.geneRemove.AS.grl = all.assembled.exon.geneRemove.grl
#strand(all.assembled.exon.geneRemove.AS.grl) = gsub("\\+", "\\:", strand(all.assembled.exon.geneRemove.AS.grl))
#strand(all.assembled.exon.geneRemove.AS.grl) = gsub("\\-", "\\+", strand(all.assembled.exon.geneRemove.AS.grl))
#strand(all.assembled.exon.geneRemove.AS.grl) = gsub("\\:", "\\-", strand(all.assembled.exon.geneRemove.AS.grl))
tmp.strand = gsub("\\+", "\\:", strand(all.assembled.exon.geneRemove.AS.grl))
tmp.strand = gsub("\\-", "\\+", tmp.strand)
tmp.strand = gsub("\\:", "\\-", tmp.strand)
strand(all.assembled.exon.geneRemove.AS.grl) = tmp.strand

all.assembled.exon.gene.AS.ol = findOverlaps(all.assembled.exon.geneRemove.AS.grl, geneModel.exon.grl)
length(unique(queryHits(all.assembled.exon.gene.AS.ol)))
all.assembled.exon.gene.AS.ol.df = data.frame(Assembled = names(all.assembled.exon.geneRemove.AS.grl)[queryHits(all.assembled.exon.gene.AS.ol)],
			      Gene = names(geneModel.exon.grl)[subjectHits(all.assembled.exon.gene.AS.ol)])
write.table(all.assembled.exon.gene.AS.ol.df, file = paste(lincRNA.path, "/Output/", filename, ".assembled.exon.gene.AS.ol.txt", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)

all.assembled.exon.gene.AS.grl = all.assembled.exon.geneRemove.grl[unique(queryHits(all.assembled.exon.gene.AS.ol))]
all.assembled.exon.gene.AS.df = as.data.frame(unlist(all.assembled.exon.gene.AS.grl), row.names = NULL)
write.table(all.assembled.exon.gene.AS.df, file = paste(lincRNA.path, "/Output/", filename, ".assembled.exon.gene.AS.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
all.assembled.exon.gene.ASremove.grl = subset(all.assembled.exon.geneRemove.grl, !(names(all.assembled.exon.geneRemove.grl) %in% names(all.assembled.exon.gene.AS.grl)))

all.assembled.exon.gene.intronic.ol = findOverlaps(all.assembled.exon.gene.ASremove.grl, geneModel.mRNA.gr)
all.assembled.exon.gene.intronic.grl = all.assembled.exon.gene.ASremove.grl[unique(queryHits(all.assembled.exon.gene.intronic.ol))]
all.assembled.exon.gene.intronicRemove.grl = subset(all.assembled.exon.gene.ASremove.grl, !(names(all.assembled.exon.gene.ASremove.grl) %in% names(all.assembled.exon.gene.intronic.grl)))

all.assembled.exon.gene.intronicRemove.AS.grl = all.assembled.exon.gene.intronicRemove.grl
#strand(all.assembled.exon.gene.intronicRemove.AS.grl) = gsub("\\+", "\\:", strand(all.assembled.exon.gene.intronicRemove.AS.grl))
#strand(all.assembled.exon.gene.intronicRemove.AS.grl) = gsub("\\-", "\\+", strand(all.assembled.exon.gene.intronicRemove.AS.grl))
#strand(all.assembled.exon.gene.intronicRemove.AS.grl) = gsub("\\:", "\\-", strand(all.assembled.exon.gene.intronicRemove.AS.grl))
tmp.strand = gsub("\\+", "\\:", strand(all.assembled.exon.gene.intronicRemove.AS.grl))
tmp.strand = gsub("\\-", "\\+", tmp.strand)
tmp.strand = gsub("\\:", "\\-", tmp.strand)
strand(all.assembled.exon.gene.intronicRemove.AS.grl) = tmp.strand

all.assembled.exon.gene.intronic.AS.ol = findOverlaps(all.assembled.exon.gene.intronicRemove.AS.grl, geneModel.mRNA.gr)
all.assembled.exon.gene.intronic.AS.grl = all.assembled.exon.gene.intronicRemove.grl[unique(queryHits(all.assembled.exon.gene.intronic.AS.ol))]

all.assembled.exon.gene.intergenic.grl = subset(all.assembled.exon.gene.intronicRemove.grl, !(names(all.assembled.exon.gene.intronicRemove.grl) %in% names(all.assembled.exon.gene.intronic.AS.grl)))
#all.assembled.exon.gene.intergenic.grl = all.assembled.exon.gene.intronicRemove.grl[-unique(queryHits(all.assembled.exon.gene.intronic.AS.ol))] 
all.assembled.exon.gene.intergenic.df = as.data.frame(unlist(all.assembled.exon.gene.intergenic.grl), row.names = NULL)
write.table(all.assembled.exon.gene.intergenic.df, file = paste(lincRNA.path, "/Output/", filename, ".assembled.exon.gene.intergenic.txt", sep = ""), quote = F, sep = "\t", row.names = F, col.names = T)

all.assembled.num = c(
		  length(all.assembled.exon.all.grl),
		  length(all.assembled.exon.grl),
		  length(geneModel.exon.grl),
		  length(all.assembled.exon.gene.grl), 
		  length(unique(subjectHits(all.assembled.exon.gene.ol))),
		  length(all.assembled.exon.gene.AS.grl),
		  length(unique(subjectHits(all.assembled.exon.gene.AS.ol))),
		  length(all.assembled.exon.gene.intronic.grl), 
		  length(all.assembled.exon.gene.intronic.AS.grl),
		  length(all.assembled.exon.gene.intergenic.grl))

names(all.assembled.num) = c("assembed", "assembled_long", "gene", "assembled_gene", "gene_assembled", "assembled_AS", "gene_AS", "assembled_intronic", "assembled_ASintronic", "assembled_intergenic")
write.table(data.frame(all.assembled.num), file = paste(lincRNA.path, "/Output/", filename, ".assembled.num.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = F)

save(all.assembled.exon.gene.grl,
     all.assembled.exon.gene.AS.grl,
     all.assembled.exon.gene.intronic.grl, 
     all.assembled.exon.gene.intronic.AS.grl,
     all.assembled.exon.gene.intergenic.grl, 
     file = paste(lincRNA.path, "/Rdata/", filename, ".assembled.annotation.grl.RData", sep = ""))
}

files = list.files(path = "../../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      myfun(fname = files[i])
}

print(paste("The program ends at: ", Sys.time()))
