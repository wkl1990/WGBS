###
#Descripition: call snp using Bis-SNP
#Environment: shell
#Author: Wang Kangli (CSU)
#Date: 2018-07-13
###

#set file
input="R1_val_1_bismark_bt2_pe.deduplicated.bam"
sortbam="R1_val_1_bismark_bt2_pe.deduplicated.sorted.bam" 
groupbam="R1_val_1_bismark_bt2_pe.deduplicated.sorted.group.bam"
methvct="R1.meth.vcf"
snpvcf="R1.snp.vcf"
snpsortvcf="R1.snp.sorted.vcf"
snpfiltervcf="R1.snp.sorted.filtered.vcf"
snpfilterlog="R1.snp.sorted.filter.summary.txt"

#sort bam file and index
samtools sort $input $sortbam
samtools index $sortbam

#add a new read group to the library
java -jar /picard-tools-1.119/AddOrReplaceReadGroups.jar I=$sortbam O=$groupbam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

#call snp (attention parameters: nonDirectional for non directional library protocol)
java -Xmx24g -jar /BisSNP/BisSNP-0.82.2.jar -R /genomes/human_UCSC_hg19/chrAll.fasta -I $groupbam -T BisulfiteGenotyper --trim_5_end_bp 0 --trim_3_end_bp 0 -vfn1 $methvct -vfn2 $snpvcf -mbq 12 -minConv 1 -toCoverage 1000 -mmq 20 --dbsnp /dbsnp/dbsnp147.GRCH37/All_20160601.vcf.gz -nt 4 # -L chr11:7000000-7100000

#sort vcf and filter fake SNPs (By default, it filter out SNPs with quality score less than 20, reads coverage more than 250, strand bias more than -0.02, quality score by depth less than 1.0, mapping quality zero reads fraction more than 0.1 and 2 SNPs within the same 10 bp window.)
perl /BisSNP/sortByRefAndCor.pl --k 1 --c 2 --tmp bissnp/tmp $snpvcf /genomes/human_UCSC_hg19/chrAll.fasta.fai > $snpsortvcf
java -Xmx24g -jar /BisSNP/BisSNP-0.82.2.jar -R genomes/human_UCSC_hg19/chrAll.fasta -T VCFpostprocess -oldVcf $snpsortvcf -newVcf $snpfiltervcf -snpVcf $snpsortvcf -o $snpfilterlog
