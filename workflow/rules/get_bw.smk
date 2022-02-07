## convert the bam file to bigwig format (signal normalized to CPM) for viewing the data on Genome Browse

rule get_bw:
    input: 
        rmBList = "STAR_align_filter/{sample}.rmBList.bam"
    output: 
        bw_cpm = "bw_cpm/{sample}.bw"
    log: "logs/get_bw/{sample}.log"
    threads: 8
    params:
        mem = '10G',
        jobName = "get_bw.{sample}"
    shell:
        ## strand of the data is opposite to dUTP method, as the UMI adapter process flipped the strandedness 
        '''
        bamCoverage -b {input} -o {output.bw_cpm} --exactScaling --normalizeUsing CPM --binSize 10 \\
        --numberOfProcessors {threads} \\
        --blackListFileName /home/dic/public/genomes/UCSC/mm10/mm10-blacklist.v2.bed &> {log}
        '''

