## ----setup, include=FALSE, echo = FALSE, cache = FALSE, warning=FALSE---------
# replace with path where you want the results be
knitr::opts_knit$set(root.dir="/home/dic/bassing_lab/users/dic/ATAC-seq/workflow/")


## ----include=FALSE,echo=FALSE,message = FALSE---------------------------------
library(ATACseqQC)
# library(BSgenome.Mmusculus.UCSC.mm10)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)


## ----message=FALSE, warning=FALSE, fig.height=6, fig.width=6------------------
# setwd("/home/dic/bassing_lab/users/dic/ATAC-seq/workflow/")
## input is bamFile
#bamfile <- list.files(path="./STAR_align_filter", pattern=glob2rx("ctrlA.rmBList.bam$"), full.names=TRUE)
bamfile <- snakemake@input[["rmBList"]]
bamfile.labels <- sub(".rmBList.bam", "", basename(bamfile))

## generate fragement size distribution
pdf(snakemake@output[["fragsize"]])
fragSizeDist(bamfile, bamfile.labels)
dev.off()

## test
# bamfile <- system.file("extdata", "GL1.bam", package="ATACseqQC", mustWork=TRUE)
# bamfile.labels <- gsub(".bam", "", basename(bamfile))
# pdf("fragsize.pdf")
# fragSizeDist(bamfile, bamfile.labels)
# dev.off()

