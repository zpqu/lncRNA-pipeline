##Make tables and plots for annotation of lincRNAs 

##Step 1, get summary table for all reconstructed transcripts

setwd("../../Output")
rm(list = ls())
num.files = list.files("../Output/", pattern = ".num.txt$")
num.all.df = as.data.frame(matrix(0, 10, length(num.files)))
for(i in 1:length(num.files)){
      num.name = gsub("\\.assembled\\.num\\.txt", "", num.files[i])
      num = read.delim(paste("../Output/", num.files[i], sep = ""), header = F)
      num.all.df[,i] = num$V2
      colnames(num.all.df)[i] = num.name
      rownames(num.all.df) = num$V1
}

write.table(num.all.df, file = "../Summary/all.assembled.num.txt", quote = F, sep = "\t", col.names = T, row.names = T)

##Step 2, get summary table for long intergenic ncRNAs
library(Rsamtools)

sam.files = num.files
sam.files = gsub("\\.assembled\\.num\\.txt", "", sam.files)
intergenic.num.df = as.data.frame(matrix(0, 8, length(num.files)))
for(i in 1:length(sam.files)){
      all.seq = as.numeric(system(paste("grep '>' ../Output/", sam.files[i], ".assembled.exon.gene.intergenic.fa | wc -l", sep = ""), intern = T))
      noBlast.seq = as.numeric(system(paste("grep '>' ../ORF/", sam.files[i], ".assembled.exon.gene.intergenic.noBlastx.fa | wc -l", sep = ""), intern = T))
      blast.seq = all.seq - noBlast.seq
      noORF.seq = as.numeric(system(paste("grep '>' ../ORF/", sam.files[i], ".assembled.exon.gene.intergenic.noBlastx.noORF.fa | wc -l", sep = ""), intern = T))
      orf.seq = noBlast.seq - noORF.seq
      classify.lincRNA = read.delim(paste("../LincRNA/", sam.files[i], ".lincRNA.neigh.annotated.txt", sep = ""), header = F)
      classify.all = paste(classify.lincRNA[,ncol(classify.lincRNA)-1], classify.lincRNA[, ncol(classify.lincRNA)], sep = "")     
      classify.all = classify.all[-grep("NA", classify.all)]
      classify.lincRNA.num = as.numeric(table(classify.all))
      intergenic.num.df[,i] = c(all.seq, blast.seq, orf.seq, noORF.seq, classify.lincRNA.num)
      colnames(intergenic.num.df)[i] = sam.files[i]
}

rownames(intergenic.num.df) = c("Intergenic", "Blastx", "ORF", "LincRNA", "AS3END", "AS5END", "Sense3END", "Sense5END")
write.table(intergenic.num.df, file = "../Summary/intergenic.assembled.num.txt", quote = F, sep = "\t", col.names = T, row.names = T)
