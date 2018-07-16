###
#Descripition: allele specific methylation
#Environment: shell
#Author: Wang Kangli (CSU)
#Date: 2018-07-16
###

#using methpipe 
#set file
bamfile="R1_val_1_bismark_bt2_pe.deduplicated.bam"
mrfile="R1_val_1_bismark_bt2_pe.deduplicated.mr"
mrsortfile="R1_val_1_bismark_bt2_pe.deduplicated.mr.sorted_start"
epireadfile="R1_val_1_bismark_bt2_pe.deduplicated.epiread"
amrfile="R1_val_1_bismark_bt2_pe.deduplicated.amr"
#format transformation 
to-mr -o $mrfile -m bismark $bamfile
LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o $mrsortfile $mrfile
methstates -c /genomes/human_UCSC_hg19/chromFa/ -o $epireadfile $mrsortfile
#Allelically methylated regions (AMRs) (attention parameters: w for size of sliding window, b to use BIC to compare models)
amrfinder -o $amrfile -c /genomes/human_UCSC_hg19/chrAll.fasta $epireadfile

#using bismark
 