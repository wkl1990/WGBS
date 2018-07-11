#Mapping, deduplication, and call methylation using bismark for paired-end reads
#transform genome reference (run only the first time, you must get the genome reference first)
bismark_genome_preparation --path_to_bowtie /bowtie2-2.2.6 --verbose /Reference/human_UCSC_hg19/
#set the file
read1="R1_val_1.fq.gz"
read2="R2_val_2.fq.gz"
mappinglogfile="mapping.log"
deduplogfile="dedup.log"
methlogfile="meth.log"
mapfile="R1_val_1_bismark_bt2_pe.bam"
dedupfile="R1_val_1_bismark_bt2_pe.deduplicated.bam"
#mapping (attention parameter: se for single-end, phred33-quals or phred64-quals, non_directional and pbat for specific library)
bismark --bowtie2 -p 4 --temp_dir bismark/temp_dir -o bismark/output_dir /Reference/human_UCSC_hg19/ -1 $read1 -2 $read2 2>$mappinglogfile
#dedup (not recommended for RRBS, attention parameter: s for single-end and p for paired-end)
deduplicate_bismark -p --bam $mapfile 2>$deduplogfile
#extraction (attention parameter: s for single-end and p for paired-end, cutoff for minimum number coverage, bedGraph to get bedGraph file, cytosine_report to get methylation report)
bismark_methylation_extractor -p --no_overlap --comprehensive --gzip --parallel 4 --bedGraph --cufoff 5 --cytosine_report -o bismark/meth_dir --genome_folder /Reference/human_UCSC_hg19/ $dedupfile 2>$methlogfile
