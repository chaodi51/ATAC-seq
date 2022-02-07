rule trim_galore:
    input: 
        get_fastq ## call the function in main Snakefile
    output:
        fq1 = "trimmed_fq/{sample}_1_trimmed.fq.gz", # always ends with fq.gz not fastq.gz
        fq1_rep = "trimmed_fq/{sample}_1.fastq.gz_trimming_report.txt", # "_1.fastq.gz" is from the original raw reads file, sometime it is "fq.gz"
        fq2 = "trimmed_fq/{sample}_2_trimmed.fq.gz",
        fq2_rep = "trimmed_fq/{sample}_2.fastq.gz_trimming_report.txt"
    log: "logs/trimmed/{sample}.log"
    threads: 10
    params:
        outdir = trimmed_fq,
        mem = '6G',
        jobName = "trim_galore.{sample}" 
    shell:
        '''
        trim_galore --cores {threads} --gzip --fastqc --paired -o {params.outdir} {input} &> {log} \
        mv trimmed_fq/{wildcards.sample}_1_val_1.fq.gz {output.fq1} &>> {log} \
        mv trimmed_fq/{wildcards.sample}_2_val_2.fq.gz {output.fq2} &>> {log} 
        '''

rule multiqc_fastqc:
    input: expand(['trimmed_fq/{sample}_1_val_1_fastqc.html','trimmed_fq/{sample}_2_val_2_fastqc.html'], sample=SAMPLES)
    output: '../results/multiqc/fastqc/multiqc_report.html'
    params:
        mem = '6G',
        jobName = "multiqc_fastqc",
        outdir = "../results/multiqc/fastqc"
    log: 'logs/multiqc/fastqc/multiqc.log'

    shell:
        '''
        rm -r {params.outdir}
        multiqc trimmed_fq/ -o {params.outdir} &> {log}
        '''