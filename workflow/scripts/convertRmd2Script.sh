
## convert .rmd file to a R script
Rscript -e "f <- 'ATACseqQC.R'; knitr::purl('ATACseqQC.Rmd', output=f)"
