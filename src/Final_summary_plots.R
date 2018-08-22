##Make tables and plots for annotation of lincRNAs 

##Step 1, get summary table for all reconstructed transcripts

rm(list = ls())
num.files = list.files("../../Output/", pattern = ".num.txt$")
num.all.df = as.data.frame(matrix(0, 10, 3))
for(i in 1:length(num.files)){
      num.name = gsub("\\.assembled\\.num\\.txt", "", num.files[i])
      num = read.delim(paste("../../Output/", num.files[i], sep = ""), header = F)
      num.all.df[,i] = num$V2
      colnames(num.all.df)[i] = num.name
      rownames(num.all.df) = num$V1
}

write.table(num.all.df, file = "../../Summary/all.assembled.num.txt", quote = F, sep = "\t", col.names = T, row.names = T)

##Step 2, get summary table for intergenic, exon_AS, intron_AS, intronic_AS
library(Rsamtools)

sam.files = num.files
sam.files = gsub("\\.assembled\\.num\\.txt", "", sam.files)
intergenic.num.df = as.data.frame(matrix(0, 8, 3))
exonAS.num.df = intronAS.num.df = intronic.num.df = as.data.frame(matrix(0, 4, 3))
for(i in 1:length(sam.files)){
      all.seq = as.numeric(system(paste("grep '>' ../../Output/", sam.files[i], ".assembled.exon.gene.intergenic.fa | wc -l", sep = ""), intern = T))
      noBlast.seq = as.numeric(system(paste("grep '>' ../../ORF/", sam.files[i], ".assembled.exon.gene.intergenic.noBlastx.fa | wc -l", sep = ""), intern = T))
      blast.seq = all.seq - noBlast.seq
      noORF.seq = as.numeric(system(paste("grep '>' ../../ORF/", sam.files[i], ".assembled.exon.gene.intergenic.noBlastx.noORF.fa | wc -l", sep = ""), intern = T))
      orf.seq = noBlast.seq - noORF.seq
      classify.lincRNA = read.delim(paste("../../LincRNA/", sam.files[i], ".lincRNA.neigh.annotated.txt", sep = ""), header = F)
      classify.all = paste(classify.lincRNA[, ncol(classify.lincRNA)-1], classify.lincRNA[, ncol(classify.lincRNA)], sep = "")     
      classify.all = classify.all[complete.cases(classify.lincRNA)]
      classify.lincRNA.num = as.numeric(table(classify.all))
      intergenic.num.df[,i] = c(all.seq, blast.seq, orf.seq, noORF.seq, classify.lincRNA.num)
      colnames(intergenic.num.df)[i] = sam.files[i]

      exonAS.all.seq = as.numeric(system(paste("grep '>' ../../Output/", sam.files[i], ".assembled.exon.gene.AS.fa | wc -l", sep = ""), intern = T))
      exonAS.noBlast.seq = as.numeric(system(paste("grep '>' ../../ORF/", sam.files[i], ".assembled.exon.gene.AS.noBlastx.fa | wc -l", sep = ""), intern = T))
      exonAS.blast.seq = exonAS.all.seq - exonAS.noBlast.seq
      exonAS.noORF.seq = as.numeric(system(paste("grep '>' ../../ORF/", sam.files[i], ".assembled.exon.gene.AS.noBlastx.noORF.fa | wc -l", sep = ""), intern = T))
      exonAS.orf.seq = exonAS.noBlast.seq - exonAS.noORF.seq
      exonAS.num.df[,i] = c(exonAS.all.seq, exonAS.blast.seq, exonAS.orf.seq, exonAS.noORF.seq)
      colnames(exonAS.num.df)[i] = sam.files[i]

      intronAS.all.seq = as.numeric(system(paste("grep '>' ../../Output/", sam.files[i], ".assembled.exon.gene.intronic.AS.fa | wc -l", sep = ""), intern = T))
      intronAS.noBlast.seq = as.numeric(system(paste("grep '>' ../../ORF/", sam.files[i], ".assembled.exon.gene.intronic.AS.noBlastx.fa | wc -l", sep = ""), intern = T))
      intronAS.blast.seq = intronAS.all.seq - intronAS.noBlast.seq
      intronAS.noORF.seq = as.numeric(system(paste("grep '>' ../../ORF/", sam.files[i], ".assembled.exon.gene.intronic.AS.noBlastx.noORF.fa | wc -l", sep = ""), intern = T))
      intronAS.orf.seq = intronAS.noBlast.seq - intronAS.noORF.seq
      intronAS.num.df[,i] = c(intronAS.all.seq, intronAS.blast.seq, intronAS.orf.seq, intronAS.noORF.seq)
      colnames(intronAS.num.df)[i] = sam.files[i]

      intronic.all.seq = as.numeric(system(paste("grep '>' ../../Output/", sam.files[i], ".assembled.exon.gene.intronic.fa | wc -l", sep = ""), intern = T))
      intronic.noBlast.seq = as.numeric(system(paste("grep '>' ../../ORF/", sam.files[i], ".assembled.exon.gene.intronic.noBlastx.fa | wc -l", sep = ""), intern = T))
      intronic.blast.seq = intronic.all.seq - intronic.noBlast.seq
      intronic.noORF.seq = as.numeric(system(paste("grep '>' ../../ORF/", sam.files[i], ".assembled.exon.gene.intronic.noBlastx.noORF.fa | wc -l", sep = ""), intern = T))
      intronic.orf.seq = intronic.noBlast.seq - intronic.noORF.seq
      intronic.num.df[,i] = c(intronic.all.seq, intronic.blast.seq, intronic.orf.seq, intronic.noORF.seq)
      colnames(intronic.num.df)[i] = sam.files[i]

}
rownames(intergenic.num.df) = c("Intergenic", "Blastx", "ORF", "LincRNA", "AS3END", "AS5END", "Sense3END", "Sense5END")
rownames(exonAS.num.df) = rownames(intronAS.num.df) = rownames(intronic.num.df) = c("All", "Blastx", "ORF", "ncRNA")

write.table(intergenic.num.df, file = "../../Summary/intergenic.assembled.num.txt", quote = F, sep = "\t", col.names = T, row.names = T)
write.table(exonAS.num.df, file = "../../Summary/exonAS.assembled.num.txt", quote = F, sep = "\t", col.names = T, row.names = T)
write.table(intronAS.num.df, file = "../../Summary/intronAS.assembled.num.txt", quote = F, sep = "\t", col.names = T, row.names = T)
write.table(intronic.num.df, file = "../../Summary/intronic.assembled.num.txt", quote = F, sep = "\t", col.names = T, row.names = T)
