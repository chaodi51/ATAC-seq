rule ATACseqQC:
    input: 
        rmBList = "STAR_align_filter/{sample}.rmBList.bam"
    output: 
        ATACseqQC_report = "ATACseqQC/{sample}_report.html"
    threads: 1
    params:
        mem = '80G',
        jobName = "ATACseqQC.{sample}"
    log: "logs/ATACseqQC/{sample}.log"
    script:
        "../scripts/ATACseqQC.Rmd"
        