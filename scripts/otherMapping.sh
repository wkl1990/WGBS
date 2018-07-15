###
#Description: mapping test using other methods
#Environment: shell
#Author: Wang Kangli (CSU)
#Date: 2018-07-13
###

#set the files
read1="R1_val_1.fq.gz"
read2="R2_val_2.fq.gz"

#bsmap 
#mapping (attention parameters: p for number of processors, v for the mismatch, g for the gap size, w for max number of equal best hits to count, n for specific library, z for base quality)
bsmap -a $read1 -b $read2 -p 8 -d /genomes/human_UCSC_hg19/chrAll.fasta -o R1.sam -v 0.1 2>bsmap.log 
bsmap -a $read1 -b $read2 -p 8 -d /genomes/human_UCSC_hg19/chrAll.fasta -v 0.1 | samtools view -bS - > R1.bam
#call methylation rations (attention parameters: u for only unique mappings, p for only properly paired mappings, z for report loci with zero methylation ratios, r to remove duplicated reads.)
methratio.py -d /genomes/human_UCSC_hg19/chrAll.fasta -o out.meth.txt -u -p -z -r R1.bam 2>methration.log

#bwameth
#genome reference index
bwameth.py index /Reference/human_UCSC_hg19/chrAll.fasta
samtools faidx /Reference/human_UCSC_hg19/chrAll.fasta
#mapping
bwameth.py --threads 8 --prefix R1 --reference Reference/human_UCSC_hg19/chrAll.fasta $read1 $read2 2>bwameth_align.log &
#check
samtools flagstat R1.bam > R1_flagstat.txt
#call methylation
bwameth.py tabulate --trim 3,3 --map-q 60 --bissnp /BisSNP/BisSNP-0.82.2.jar --prefix R1 -t 12 --reference Reference/human_UCSC_hg19/chrAll.fasta R1.bam

#WALT
#makedb reference
makedb -c /genomes/human_UCSC_hg19/chrAll.fasta -o /Reference/human_UCSC_hg19/chrAll.fasta.dbindex
#mapping
walt -t 8 -i /Reference/human_UCSC_hg19/chrAll.fasta.dbindex -1 $read1 -2 $read2 -o R1.mr
#merging libraries and removing duplicates
LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o R1.mr.sorted_start R1.mr
duplicate-remover -S R1.dremove_stat.txt -o R1.mr.dremove R1.mr.sorted_start
#Estimating bisulfite conversion rate
bsrate -c /Reference/human_UCSC_hg19/ -o R1.bsrate R1.mr.dremove
#computing single-site methylation levels using methcounts
methcounts -c /Reference/human_UCSC_hg19/ -o R1.dremove.meth R1.mr.dremove
#extracting and merging symmetric CpG methylation levels
symmetric-cpgs -o R1_CpG.meth R1.dremove.meth
#computation of methylation level statistics
levels -o R1.levels R1.dremove.meth
