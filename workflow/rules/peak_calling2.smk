rule macs2_narrow2:
    input: 
        merged_bam = "merged_bams/{group}.bam"
    output: 
        narrowPeak = "macs2_merged_bam/narrowPeak/{group}_peaks.narrowPeak"
    threads: 1
    params:
        mem = '20G',
        jobName = "macs2_narrow2.{group}",
        narrow_out = "macs2_merged_bam/narrowPeak"
    log: "logs/macs2_narrow2/{group}.log"
    shell:
        '''
        mkdir -p macs2_merged_bam/narrowPeak 
        macs2 callpeak -f BAMPE -g mm --keep-dup all -n {wildcards.group} -t {input.merged_bam} --outdir {params.narrow_out} &> {log}    
        '''

rule macs2_broad2:
    input: 
        merged_bam = "merged_bams/{group}.bam"
    output: 
        broadPeak = "macs2_merged_bam/broadPeak/{group}_peaks.broadPeak" 
    threads: 1
    params:
        mem = '20G',
        jobName = "macs2_broad2.{group}",
        broad_out = "macs2_merged_bam/broadPeak"
    log: "logs/macs2_broad2/{group}.log"
    shell:
        '''
        mkdir -p macs2_merged_bam/broadPeak
        macs2 callpeak --broad -f BAMPE -g mm --keep-dup all -n {wildcards.group} -t {input.merged_bam} --outdir {params.broad_out} &>> {log}    
        '''