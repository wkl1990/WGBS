###
#Description: Illumina paired-end reads quality control
#Environment: shell
#Author: Wang Kangli (CSU)
#Date: 2018-07-13
###

#set the files 
read1="R1.fastq.gz"
read2="R2.fastq.gz"
fastqclogfile="fastqc.log"
trimlogfile="trimgalore.log"

#fastqc
fastqc --outdir fastqc $read1 $read2 > $fastqclogfile

#trim (attention parameters: phred33 or phred64, adapter)
trim_galore --stringency 5 --clip_R1 5 --clip_R2 5 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --fastqc --fastqc_args "--outdir fastqc" -o trimed --paired $read1 $read2 > $trimlogfile
