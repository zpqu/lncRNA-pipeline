###fun


rm(list = ls())
myfun = function(file){
      gtf = read.delim(paste("../GTF/", file, ".gtf", sep = ""), header = F)
      gtf.inter = read.delim(paste("../LincRNA/", file, ".lincRNA.neigh.annotated.txt", sep = ""), header = F)
      test.id = gtf$V9
      test.id = gsub(".+transcript\\_id\\s", "", test.id)
      test.id = gsub("\\;\\s.+", "", test.id)

      gtf.inter.gtf = subset(gtf, test.id %in% gtf.inter$V6)
      write.table(gtf.inter.gtf, file = paste("../LincRNA/", file, ".lincRNA.gtf", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)
}

setwd("../../LincRNA")
files = list.files(path = "../GTF/", pattern = "gtf$")
for(i in 1:length(files)){
      fname = gsub("\\.gtf","",files[i])
      myfun(file = fname)
}
