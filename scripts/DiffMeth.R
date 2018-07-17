###
#Description: differential methylation analysis using methylKit, Comb-p, metilene
#Environment: shell
#Author: Wang Kangli (CSU)
#Date: 2018-07-16
###

#load packages
library(methylKit)
library(xlsx)
library(genomation)
library(GenomicRanges)
library(qvalue)

#Set the methylation call files
file.list <- list("S1_val_1_bismark_bt2_pe_CpG_report_merge.merged_CpG_evidence.cov","S2_val_1_bismark_bt2_pe_CpG_report_merge.merged_CpG_evidence.cov","S3_val_1_bismark_bt2_pe_CpG_report_merge.merged_CpG_evidence.cov","S4_val_1_bismark_bt2_pe_CpG_report_merge.merged_CpG_evidence.cov")
#Read the sample information (sampleID, treatment, batch and other information could be added)
sample_trait <- read.xlsx("sample_info.xlsx",sheetName = "Sheet1")
batch <- sample_trait$batch 
sampleID <- as.list(as.character(sample_trait$sampleID)) #Attention to transform to list to match the file.list.
treatment <- as.numeric(sample_trait$treatment) #Attention to transform to numeric
#Read methylation files
meth_cov <- methRead(file.list,sample.id=sampleID,assembly="hg19",pipeline="bismarkCoverage",treatment=treatment,header=FALSE,context="CpG")

#Descriptive statistics on samples (set i to number for each samples)
getMethylationStats(meth_cov[[i]],plot=TRUE,both.strands=FALSE)
getCoverageStats(meth_cov[[i]],plot=TRUE,both.strands=FALSE)

#Filtering samples based on read coverage (optional)
filtered.meth_cov <- filterByCoverage(meth_cov,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)

#Comparative analysis
#Merging samples
meth <- unite(meth_cov, destrand=FALSE)
head(meth)
#Sample correlation plot
getCorrelation(meth,plot=TRUE)
#Clustering samples
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
#PCA plot
PCASamples(meth, screeplot=TRUE)

#Batch effects (combat method is recommended, removeComp is optional)
#removeComp method
sampleAnnotation <- sample_trait[,c("treatment","gender","batch")]
as <- assocComp(mBase=meth,sampleAnnotation)
newMeth <- removeComp(meth,comp=1)
#combat method
mat <- percMethylation(meth,rowids=TRUE)
mat_combat <- ComBat(mat,batch)
newMeth <- reconstruct(mat_combat,meth)
length(grep("chr[0-9]",rownames(mat)))

#Finding differentially methylated bases (calculateDiffMeth with Chisq-test or the F-test is used, other methods such as calculateDiffMethDSS and other appropriate statistical test also can be used)
#calculateDiffMeth (attention parameters: overdispersion to correct for overdispersion, test to choose the statistical test, mc.cores to use multiple cores, covariates to account for covariates)
covariates <- sample_trait$gender
myDiff <- calculateDiffMeth(meth,covariates=covariates,overdispersion="MN",test="Chisq",mc.cores=1)
#get hyper, hypo and all methylated bases
myDiff25p.hyper <- getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo <- getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
myDiff25p <- getMethylDiff(myDiff,difference=25,qvalue=0.01)
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)


#Annotating differentially methylated bases
#read the gene BED file (customized annotation is recommoned)
gene.obj <- readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))
#annotate differentially methylated CpGs with promoter/exon/intron using annotation data
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
#read the shores and flanking regions and name the flanks as shores and CpG islands as CpGi
cpg.obj <- readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit"), feature.flank.name=c("CpGi","shores"))
#convert methylDiff object to GRanges and annotate
diffCpGann <- annotateWithFeatureFlank(as(myDiff25p,"GRanges"), cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")
#get the distance to TSS and nearest gene name 
diffAnn <- annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
head(getAssociationWithTSS(diffAnn))
#get percentage/number of differentially methylated regions that overlap with intron/exon/promoters
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
#plot the percentage of differentially methylated bases overlapping with exon/intron/promoters
plotTargetAnnotation(diffAnn,precedence=TRUE, main="differential methylation annotation")
#plot the CpG island annotation
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"), main="differential methylation annotation")
#get percentage of intron/exon/promoters that overlap with differentially methylated bases.
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

#customized annotation (annotate differentially methylated bases using genes based on the physical position and chromatin interaction)
#read reference
refseq_adj.obj=readTranscriptFeatures("/reference/hg19_refGene.bed", remove.unusual = FALSE, up.flank = 2000, down.flank = 0, unique.prom = FALSE)
refseq_adj_diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),refseq_adj.obj)
#annotate bases within promoter/intron/exon
promoter_overlap <- findOverlaps(as(myDiff25p,"GRanges"),refseq_adj.obj$promoters)
promoter_overlap_gene <- refseq_adj.obj$promoters$name[promoter_overlap@to]
promoter_overlap_all <- cbind(as.character(myDiff25p[promoter_overlap@from,1]),promoter_overlap_gene)
colnames(promoter_overlap_all) <- c("CpG","refseq")
promoter_overlap_all <- as.data.frame(promoter_overlap_all)
promoter_overlap_all$type <- "promoter"
#exon
exon_overlap <- findOverlaps(as(myDiff25p,"GRanges"),refseq_adj.obj$exons)
exon_overlap_gene <- refseq_adj.obj$exons$name[exon_overlap@to]
exon_overlap_all <- cbind(as.character(myDiff25p[exon_overlap@from,1]),exon_overlap_gene)
colnames(exon_overlap_all) <- c("CpG","refseq")
exon_overlap_all <- as.data.frame(exon_overlap_all)
exon_overlap_all$type <- "exon"
#intron
intron_overlap <- findOverlaps(as(myDiff25p,"GRanges"),refseq_adj.obj$introns)
intron_overlap_gene <- refseq_adj.obj$introns$name[intron_overlap@to]
intron_overlap_all <- cbind(as.character(myDiff25p[intron_overlap@from,1]),intron_overlap_gene)
colnames(intron_overlap_all) <- c("CpG","refseq")
intron_overlap_all <- as.data.frame(intron_overlap_all)
intron_overlap_all$type <- "intron"
#total
allDMP_overlap_gene <- rbind(promoter_overlap_all,exon_overlap_all,intron_overlap_all)
#annotate bases within 10kb flanking of genes
DMP_unannotated_gr <- as(myDiff25p,"GRanges")[-which(myDiff25p[,1] %in% allDMP_overlap_gene$CpG),]
DMP_unannotated <- myDiff25p[-which(myDiff25p[,1] %in% allDMP_overlap_gene$CpG),]
refseq_flank.obj=readFeatureFlank("/reference/hg19_refGene.bed",remove.unusual=FALSE,flank=10000,clean=FALSE,feature.flank.name=c("genebody","within"))
genebody_overlap <- findOverlaps(DMP_unannotated_gr,refseq_flank.obj$genebody)
within10kb_overlap <- findOverlaps(DMP_unannotated_gr,refseq_flank.obj$within)
aa <- function(x){res=c() for(i in 1:length(x)){ if(x[i]<=69723){res=c(res,x[i])} else{res=c(res,x[i]-69723)}} return(res)}
within10kb_overlap_gene <- refseq_adj.obj$promoters$name[aa(within10kb_overlap@to)]
within10kb_overlap_all <- cbind(as.character(DMP_unannotated[within10kb_overlap@from,1]),within10kb_overlap_gene)
colnames(within10kb_overlap_all) <- c("CpG","refseq")
within10kb_overlap_all <- as.data.frame(within10kb_overlap_all)
within10kb_overlap_all$type <- "within10kb"
#annotate bases based on chromatin interaction (now is not available)

#statistics and plot
#gene features
diffAnn_myDiff_refseq_adj=annotateWithGeneParts(as(myDiff25p,"GRanges"),refseq_adj.obj)
head(getAssociationWithTSS(diffAnn_myDiff_refseq_adj))
getTargetAnnotationStats(diffAnn_myDiff_refseq_adj,percentage=TRUE,precedence=TRUE)
getFeatsWithTargetsStats(diffAnn_myDiff_refseq_adj,percentage=TRUE)
plotTargetAnnotation(diffAnn_myDiff_refseq_adj,precedence=TRUE, main="differential methylation annotation")
#CGI
cpg.obj=readFeatureFlank("/reference/CGI_Ext_hg19.bed",feature.flank.name=c("CpGi","shores"))
diffCpGann=annotateWithFeatureFlank(mas(myDiff25p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),main="differential methylation annotation")
#Volcano plot
library(ggplot2)
myDiff$meth.diff <- myDiff$FC
r03 = ggplot(myDiff,aes(meth.diff,-1*log10(pval)))
gene_anno_total <- getAssociationWithTSS(diffAnn_myDiff)
myDiff$dist.to.feature <- abs(gene_anno_total$dist.to.feature)
myDiff$significant <- "nosig"
myDiff$significant1 <- myDiff$significant
myDiff$significant[intersect(which(myDiff$pval<0.05),which(abs(myDiff$FC) > 10))] <- "down"
myDiff$significant[intersect(which(myDiff$pval<0.05),which(myDiff$FC>10))] <- "up"
myDiff$significant1[intersect(which(myDiff$significant=="down"),which(myDiff$dist.to.feature<=2000))] <- "down_withinTSS"
myDiff$significant1[intersect(which(myDiff$significant=="up"),which(myDiff$dist.to.feature<=2000))] <- "up_withinTSS"
myDiff$significant1 <- factor(myDiff$significant1,levels=c("up_withinTSS","down_withinTSS","up","down","nosig"))
r03 + geom_point(aes(color=myDiff$significant1))

#Differentially methylated regions(DMR) analysis using methylKit
#Regional analysis (tiling windows analysis compare to base-pair resolution analysis)
tiles <- tileMethylCounts(meth,win.size=1000,step.size=1000,cov.bases=2,mc.cores=8,save.db=FALSE)
head(tiles)
tiles_mat <- percMethylation(tiles,rowids=TRUE)
#Finding differentially methylated tiles and annotation just like differentially methylated bases

#DMR analysis using Comb-p
#prepare the bed file from methylKit DMB result
myDiff_seqnames <- as.character(apply(myDiff,1,function(x) strsplit(x,"\\.")[[1]][1]))
myDiff_start <- as.numeric(apply(myDiff,1,function(x) strsplit(x,"\\.")[[1]][2]))
myDiff_end <- as.numeric(apply(myDiff,1,function(x) strsplit(x,"\\.")[[1]][3]))
myDiff_bed4 <- cbind(myDiff_seqnames,myDiff_start,myDiff_end,myDiff$pval)
write.table(myDiff_bed4,file="myDiff.bed4",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
#run in shell (attention parameters: seed for a started p-value of a region, dist for extend region if find another p-value within this dist, region-filter-p for post-filter reported regions)
comb-p pipeline -c 4 --dist 500 --step 50 --seed 0.05 --prefix output myDiff.bed4

#DMR analysis using metilene 
#prepare the file from methylKit methylation level
mat_dataframe <- as.data.frame(mat)
mat_dataframe <- mat_dataframe/100
mat_dataframe$ID <- rownames(mat_dataframe)
mat_dataframe$chr <- as.character(lapply(mat_dataframe[,"ID"],function(x) strsplit(x,"\\.")[[1]][1]))
mat_dataframe$pos <- as.character(lapply(mat_dataframe[,"ID"],function(x) strsplit(x,"\\.")[[1]][2]))
mat_metilene <- mat_dataframe[,c("chr","pos",1:4)]
colnames(mat_metilene) <- c("chr","pos","g1_S1","g1_S2","g2_S3","g2_S4")
write.table(mat_metilene,file="mat_metilene.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
#run in shell (attention parameters: M for maximum distance, m for minimum cpgs, d for minimum mean methylation difference)
metilene -a g1 -b g2 mat_metilene.txt > DMR_metilene.txt
