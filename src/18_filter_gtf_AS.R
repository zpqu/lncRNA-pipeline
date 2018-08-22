###function for retrieve gtf and fa file of lncRNAs


rm(list = ls())
myfun = function(file){
      ###load gtf file
      gtf = read.delim(paste("../GTF/", file, ".gtf", sep = ""), header = F)

      test.id = gtf$V9
      test.id = gsub(".+transcript\\_id\\s", "", test.id)
      test.id = gsub("\\;\\s.+", "", test.id)
      
      ###get exonAS
      gtf.exonAS.file = paste("../LincRNA/", file, ".exonAS.neigh.txt", sep = "")
      if(file.info(gtf.exonAS.file)$size < 2){
        gtf.exonAS = NULL
      }else{
	gtf.exonAS = read.delim(gtf.exonAS.file, header = F)
      }

      gtf.exonAS.gtf = subset(gtf, test.id %in% gtf.exonAS$V6)
      write.table(gtf.exonAS.gtf, file = paste("../LincRNA/", file, ".exonAS.gtf", sep = ""), 
        quote = F, sep = "\t", col.names = F, row.names = F)
      
      ###get intronAS
      gtf = read.delim(paste("../GTF/", file, ".gtf", sep = ""), header = F)

      gtf.intronAS.file = paste("../LincRNA/", file, ".intronAS.neigh.txt", sep = "")
      if(file.info(gtf.intronAS.file)$size < 2){
        gtf.intronAS = NULL
      }else{
	gtf.intronAS = read.delim(gtf.intronAS.file, header = F)
      }

      gtf.intronAS.gtf = subset(gtf, test.id %in% gtf.intronAS$V6)
      write.table(gtf.intronAS.gtf, file = paste("../LincRNA/", file, ".intronAS.gtf", sep = ""), 
        quote = F, sep = "\t", col.names = F, row.names = F)

      ###for intronic
      gtf = read.delim(paste("../GTF/", file, ".gtf", sep = ""), header = F)

      gtf.intronic.file = paste("../LincRNA/", file, ".intronic.neigh.txt", sep = "")
      if(file.info(gtf.intronic.file)$size < 2){
        gtf.intronic = NULL
      }else{
	gtf.intronic = read.delim(gtf.intronic.file, header = F)
      }

      gtf.intronic.gtf = subset(gtf, test.id %in% gtf.intronic$V6)
      write.table(gtf.intronic.gtf, file = paste("../LincRNA/", file, ".intronic.gtf", sep = ""), 
        quote = F, sep = "\t", col.names = F, row.names = F)
}

setwd("../../LincRNA")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}
