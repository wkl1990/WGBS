###
#Description: Mapping, deduplication, and call methylation using bismark for paired-end reads
#Environment: shell
#Author: Wang Kangli (CSU)
#Date: 2018-07-13
###

#transform genome reference (run only the first time, you must get the genome reference first)
bismark_genome_preparation --path_to_bowtie /bowtie2-2.2.6 --verbose /Reference/human_UCSC_hg19/

#set the files
read1="R1_val_1.fq.gz"
read2="R2_val_2.fq.gz"
mappinglogfile="mapping.log"
deduplogfile="dedup.log"
methlogfile="meth.log"
mapfile="R1_val_1_bismark_bt2_pe.bam"
dedupfile="R1_val_1_bismark_bt2_pe.deduplicated.bam"

#mapping (attention parameters: se for single-end, phred33-quals or phred64-quals, non_directional and pbat for specific library)
bismark --bowtie2 -p 4 --temp_dir bismark/temp_dir -o bismark/output_dir /Reference/human_UCSC_hg19/ -1 $read1 -2 $read2 2>$mappinglogfile

#dedup (not recommended for RRBS, attention parameters: s for single-end and p for paired-end)
deduplicate_bismark -p --bam $mapfile 2>$deduplogfile

#extraction (attention parameters: s for single-end and p for paired-end, cutoff for minimum number coverage, bedGraph to get bedGraph file, cytosine_report to get methylation report)
bismark_methylation_extractor -p --no_overlap --comprehensive --gzip --parallel 4 --bedGraph --cufoff 5 --cytosine_report -o bismark/meth_dir --genome_folder /Reference/human_UCSC_hg19/ $dedupfile 2>$methlogfile

#get the mapping report and summary plots
#get files from previous codes
alignment_report="R1_val_1_bismark_bt2_PE_report.txt"
dedup_report="R1_val_1_bismark_bt2_pe.deduplication_report.txt"
splitting_report="R1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt"
mbias_report="R1_val_1_bismark_bt2_pe.deduplicated.M-bias.txt"
#bismark2report (If not specified files, bismark2report attempts to find report files in the current directory.)
bismark2report --alignment_report $alignment_report --dedup_report $dedup_report --splitting_report $splitting_report --mbias_report $mbias_report
#bismark2summary (If no BAM files are specified explicitly the current working directory is scanned for Bismark alignment files and their associated reports.)
bismark2summary $mapfile $dedupfile
