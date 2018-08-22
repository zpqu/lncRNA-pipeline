####This bash script is used to filter lincRNAs from all assembeled/reconstructed transcripts
####Put the gtf format of assembled/reconstructed transcripts files in to GTF/ folder
####Output will be in /LincRNA folder
####Date: 01/04/2014
####Author: Zhipeng

#!/bin/bash

echo -ne "The whole program starts at: "
date

##Step1: annotate assembed transcripts regarding to protein-coding genes
echo -ne "Now is processing Step 1 ... at "
date
cd ./src
Rscript 1_annotation_assemble_long.R

##Step2: retrieve fasta sequence for all intergenic transcripts
echo -ne "Now is processing Step 2 ... at "
date
Rscript 2_get_fasta_intergenic.R

##Step3: run blastx for all intergenic transcripts against plant UNIPROT_SPROT database
echo -ne "Now is processing Step 3 ... at "
date
#./3_run_blastx_strand.sh > ../log/3_run_blast_strand.log

##Step4: remove intergenic transcripts with hit in UNIPROT_SPROT database
echo -ne "Now is processing Step 4 ... at "
date
Rscript 4_filter_blastx.R

##Step5: find ORFs in intergenic transcripts
echo -ne "Now is processing Step 5 ... at "
date
./5_run_ORF_strand.sh >../log/5_run_ORF_strand.log

##Step6: remove intergenic transcripts with long ORFs (>100aa or > 50aa at ends)
echo -ne "Now is processing Step 6 ... at "
date
Rscript 6_filter_orf.R

##Step7: retrieve fasta sequence of lincRNAs
echo -ne "Now is processing Step 7 ... at "
date
Rscript 7_get_lincRNA.R

##Step8: find neighbour genes of lincRNAs
echo -ne "Now is processing Step 8 ... at "
date
./8_run_lincRNA_neigh.sh

##Step9: get gtf files of lincRNAs
echo -ne "Now is processing Step 9 ... at "
date
Rscript 9_filter_intergenic_gtf.R

##Step10: make summary table
echo -ne "Now is processing Step 10 ... at "
date
Rscript 10_make_summary_table.R

##Step11: get fasta file of exonAS, intronAS and intronic
echo -ne "Now is processing Step 11 ... at "
date
Rscript 11_get_fasta_AS.R

##Step12: run blastx for exonAS, intronAS and intronic
echo -ne "Now is processing Step 12 ... at "
date
#./12_run_blastx_AS_strand.sh > ../log/12_run_blastx_AS_strand.log

##Step13: filter trancripts with blastx hit(s)
echo -ne "Now is processing Step 13 ... at "
date
Rscript 13_filter_blastx_AS.R

##Step14: run getorf to find long ORFs
echo -ne "Now is processing Step 14 ... at "
date
./14_run_ORF_AS_strand.sh > ../log/14_run_ORF_AS_strand.log

##Step15: filter trancripts with long ORFs
echo -ne "Now is processing Step 15 ... at "
date
Rscript 15_filter_orf_AS.R 

##Step16: get exonAS, intronAS and intronic
echo -ne "Now is processing Step 16 ... at "
date
Rscript 16_get_AS.R

##Step17: annotate exonAS, intronAS and intronic
echo -ne "Now is processing Step 17 ... at "
date
./17_run_AS_neigh.sh

##Step18: get gtf files of exonAS, intronAS and intronic
echo -ne "Now is processing Step 18 ... at "
date
Rscript 18_filter_gtf_AS.R

##Step19: get gtf/fa files of all lncRNAs
echo -ne "Now is processing Step 19 ... at "
date
./19_merge_lncRNA.sh

##Step20: get final summary table (This script will regenate the summary table in step 10)
echo -ne "Now is processing Step 20 ... at "
date
Rscript Final_summary_plots.R

echo -ne "The whole program ends at: "
date
