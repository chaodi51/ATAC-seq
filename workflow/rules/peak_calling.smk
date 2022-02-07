rule macs2_narrow:
    input: 
        rmBList = "STAR_align_filter/{sample}.rmBList.bam"
    output: 
        narrowPeak = "macs2/narrowPeak/{sample}_peaks.narrowPeak"
    threads: 1
    params:
        mem = '10G',
        jobName = "macs2_narrow.{sample}",
        narrow_out = "macs2/narrowPeak"
    log: "logs/macs2_narrow/{sample}.log"
    shell:
        '''
        mkdir -p macs2/narrowPeak 
        macs2 callpeak -B -f BAMPE -g mm --keep-dup all -n {wildcards.sample} -t {input.rmBList} --outdir {params.narrow_out} &> {log}    
        '''

rule macs2_broad:
    input: 
        rmBList = "STAR_align_filter/{sample}.rmBList.bam"
    output: 
        broadPeak = "macs2/broadPeak/{sample}_peaks.broadPeak" 
    threads: 1
    params:
        mem = '10G',
        jobName = "macs2_broad.{sample}",
        broad_out = "macs2/broadPeak"
    log: "logs/macs2_broad/{sample}.log"
    shell:
        '''
        mkdir -p macs2/broadPeak
        macs2 callpeak -B --broad -f BAMPE -g mm --keep-dup all -n {wildcards.sample} -t {input.rmBList} --outdir {params.broad_out} &>> {log}    
        '''