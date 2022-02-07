# This script calculates RPF expression on CDS
library(Rsubread)
library(dplyr)
library(mgsub)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## run all bams together
bamfiles <- snakemake@input[["bamfiles"]]
# bamfiles <- paste0("./STAR_align_filter/", as.vector(samples$sample),".rmBList.bam")

## run one bam file
# bamfiles <- snakemake@input[["bamfile"]]
RPFcounts <- featureCounts(files=bamfiles, 
    annot.ext=snakemake@input[['merged_peaks_saf']],
    annot.inbuilt = "mm10",
    isGTFAnnotationFile=FALSE, 
    minOverlap = 1,
    fracOverlap=snakemake@params[["fracOverlap"]], 
    fracOverlapFeature=snakemake@params[["fracOverlapFeature"]],
    allowMultiOverlap=TRUE,
    strandSpecific=snakemake@params[["strand"]], 
    countMultiMappingReads=FALSE, 
    juncCounts=TRUE, 
    isPairedEnd=TRUE,
    nthreads=snakemake@threads[[1]])

write.table(RPFcounts$counts, file=snakemake@output[["featureCount"]], sep="\t", quote=F, row.names = TRUE, col.names = NA)


