# Whole genome bisulfite sequencing (WGBS) analysis
## Description
The codes contain the basic processes for WGBS analysis, which contain QC, mapping, and differential analysis. Automatic pipeline and large samples should be considered in the future.
## Softwares and packages
* Quality control: [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trim_Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* The major mapping software used: [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* Other mapping softwares: [BSMAP](https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/bsmap-for-methylation), [Bwameth](https://github.com/brentp/bwa-meth), [WALT](https://github.com/smithlabcode/walt).
Some comparison between Bismark and BSMAP in Chinese could be found [here](http://www.biotrainee.com/thread-1510-1-1.html).
* Differentially methylated analysis: [methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html), [MethPipe](http://smithlabresearch.org/software/methpipe/), [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/), [Comb-p](https://github.com/brentp/combined-pvalues)
* Allele-specific methylation: [MethPipe](http://smithlabresearch.org/software/methpipe/)
* SNP calling: [Bis-SNP](http://people.csail.mit.edu/dnaase/bissnp2011/)
* Annotation: [methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html) and customized strategy.
## Pipeline
![WGBS pipeline](http://f.cl.ly/items/343k422F1s1Y443w3g3L/WGBS.jpg)
