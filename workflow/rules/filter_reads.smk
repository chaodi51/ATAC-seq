## filter chrM and PCR duplicates
rule rmchrM:
    input: 
        bam = "STAR_align/{sample}.bam"
    output: 
        rmChrM = "STAR_align_filter/{sample}.rmChrM.bam"
    threads: 8
    params:
        mem = '50G',
        jobName = "rmchrM.{sample}"
    log: "logs/rmchrM/{sample}.log"
    shell:
        '''
        samtools view --threads {threads} -h {input.bam} | grep -v chrM | \
        samtools sort --threads {threads} -O bam -o {output.rmChrM} &> {log}
        '''
rule rmdup:
    input: 
        rmChrM = "STAR_align_filter/{sample}.rmChrM.bam"
    output: 
        rmdup = "STAR_align_filter/{sample}.rmdup.bam",
        dups = "STAR_align_filter/{sample}.dups.txt"
    threads: 8
    params:
        mem = '50G',
        jobName = "rmdup.{sample}"
    log: "logs/rmdup/{sample}.log"    
    shell:
        ''' 
        java -XX:ParallelGCThreads={threads} -jar ~/public/tools/picard-2.25.1-1/picard.jar MarkDuplicates --INPUT \
        {input.rmChrM} --OUTPUT {output.rmdup} --METRICS_FILE {output.dups} --REMOVE_DUPLICATES true &> {log}
        '''

rule filter_blacklist:
    input: 
        rmdup = "STAR_align_filter/{sample}.rmdup.bam"
    output: 
        rmBList = "STAR_align_filter/{sample}.rmBList.bam",
        index = "STAR_align_filter/{sample}.rmBList.bam.bai" 
    threads: 8
    params:
        mem = '50G',
        jobName = "filter_blacklist.{sample}",
        blacklist = config['blacklist']
    log: "logs/filter_blacklist/{sample}.log"
    shell:
        '''
        bedtools intersect -v -abam {input.rmdup} -b {params.blacklist} > {output.rmBList}

        samtools index {output.rmBList} {output.index}
        '''