####This bash script is used to filter lincRNAs from all assembeled/reconstructed transcripts
####Put the gtf format of assembled/reconstructed transcripts files in to GTF/ folder
####Output will be in /LincRNA folder
####Date: 01/04/2014
####Author: Zhipeng

cd ./src
##Step1: annotate assembed transcripts regarding to protein-coding genes
Rscript 1_annotation_assemble_long.R

##Step2: retrieve fasta sequence for all intergenic transcripts
Rscript 2_get_fasta_intergenic.R

##Step3: run blastx for all intergenic transcripts against plant UNIPROT_SPROT database
./3_run_blastx.sh > ../log/3_run_blast.log

##Step4: remove intergenic transcripts with hit in UNIPROT_SPROT database
Rscript 4_filter_blastx.R

##Step5: find ORFs in intergenic transcripts
./5_run_ORF.sh > ../log/5_run_ORF.log

##Step6: remove intergenic transcripts with long ORFs (>100aa or > 50aa at ends)
Rscript 6_filter_orf.R

##Step7: retrieve fasta sequence of lincRNAs
Rscript 7_get_lincRNA.R

##Step8: find neighbour genes of lincRNAs
./8_run_lincRNA_neigh.sh

##Step9: get gtf files of lincRNAs
Rscript 9_filter_intergenic_gtf.R

##Step10: make summary table
Rscript 10_make_summary_table.R

##Step11: get fasta file of exonAS, intronAS and intronic
Rscript 11_get_fasta_AS.R

##Step12: run blastx for exonAS, intronAS and intronic
./12_run_blastx_AS.sh > ../log/12_run_blastx_AS.log

##Step13: filter trancripts with blastx hit(s)
Rscript 13_filter_blastx_AS.R

##Step14: run getorf to find long ORFs
./14_run_ORF_AS.sh > ../log/14_run_ORF_AS.log

##Step15: filter trancripts with long ORFs
Rscript 15_filter_orf_AS.R 

##Step16: get exonAS, intronAS and intronic
Rscript 16_get_AS.R

##Step17: annotate exonAS, intronAS and intronic
./17_run_AS_neigh.sh

##Step18: get gtf files of exonAS, intronAS and intronic
Rscript 18_filter_gtf_AS.R

##Step19: merge gtf/fa files of all lncRNAs
./19_merge_lncRNA.sh

##Step20: get final summary table (This script will regenate the summary table in step 10)
Rscript Final_summary_plots.R