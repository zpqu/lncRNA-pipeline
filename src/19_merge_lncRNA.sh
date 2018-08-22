## This script is used to run blastx for intergenic transcripts against Swiss-port
## Date: 17/02/2014
## Author: Zhipeng

echo -ne "Now is merging gtf and fasta files for lncRNAs ..."
date

cd ../../GTF/
FILES=$(ls *.gtf)

cd ../LincRNA/
for f in $FILES
do
    echo -ne "Processing $f ..."
    date
    filename=${f%.*}
    sed 's/XLOC/LINC.XLOC/g' ${filename}.lincRNA.gtf >${filename}.lincRNA.rename.gtf
    sed -i 's/CUFF/LINC.CUFF/g' ${filename}.lincRNA.rename.gtf
    sed -i 's/gene/LINC/g' ${filename}.lincRNA.rename.gtf
    sed 's/CUFF/LINC.CUFF/g' ${filename}.assembled.exon.gene.lincRNA.fa >${filename}.assembled.exon.gene.lincRNA.rename.fa
    sed -i 's/gene/LINC/g' ${filename}.assembled.exon.gene.lincRNA.rename.fa

    sed 's/XLOC/exonAS.XLOC/g' ${filename}.exonAS.gtf >${filename}.exonAS.rename.gtf
    sed -i 's/CUFF/exonAS.CUFF/g' ${filename}.exonAS.rename.gtf
    sed -i 's/gene/exonAS/g' ${filename}.exonAS.rename.gtf
    sed 's/CUFF/exonAS.CUFF/g' ${filename}.assembled.exon.gene.exonAS.fa >${filename}.assembled.exon.gene.exonAS.rename.fa
    sed -i 's/gene/exonAS/g' ${filename}.assembled.exon.gene.exonAS.rename.fa

    sed 's/XLOC/intronAS.XLOC/g' ${filename}.intronAS.gtf >${filename}.intronAS.rename.gtf
    sed -i 's/CUFF/intronAS.CUFF/g' ${filename}.intronAS.rename.gtf
    sed -i 's/gene/intronAS/g' ${filename}.intronAS.rename.gtf
    sed 's/CUFF/intronAS.CUFF/g' ${filename}.assembled.exon.gene.intronAS.fa >${filename}.assembled.exon.gene.intronAS.rename.fa
    sed -i 's/gene/intronAS/g' ${filename}.assembled.exon.gene.intronAS.rename.fa

    sed 's/XLOC/INTRONIC.XLOC/g' ${filename}.intronic.gtf >${filename}.intronic.rename.gtf
    sed -i 's/CUFF/INTRONIC.CUFF/g' ${filename}.intronic.rename.gtf
    sed -i 's/gene/INTRONIC/g' ${filename}.intronic.rename.gtf
    sed 's/CUFF/INTRONIC.CUFF/g' ${filename}.assembled.exon.gene.intronic.fa >${filename}.assembled.exon.gene.intronic.rename.fa
    sed -i 's/gene/INTRONIC/g' ${filename}.assembled.exon.gene.intronic.rename.fa

    cat ${filename}.lincRNA.rename.gtf ${filename}.exonAS.rename.gtf ${filename}.intronAS.rename.gtf ${filename}.intronic.rename.gtf >${filename}.lncRNA.rename.gtf
    cat ${filename}.assembled.exon.gene.lincRNA.rename.fa ${filename}.assembled.exon.gene.exonAS.rename.fa ${filename}.assembled.exon.gene.intronAS.rename.fa ${filename}.assembled.exon.gene.intronic.rename.fa >${filename}.assembled.exon.gene.lncRNA.rename.fa    
    
done

echo -ne "The process of mergeing lncRNAs ends at: "
date
