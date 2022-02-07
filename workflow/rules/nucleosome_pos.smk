rule nucleosome_pos:
    input:
        merged_bam = "merged_bams/{group}.bam",
        bam_index = "merged_bams/{group}.bam.bai"
    output:
        pos = "nucleosome_pos/{group}_peaks.gappedPeak"
    threads: 1
    params:
        mem = '100G',
        jobName = "nucleosome_pos.{group}",
        genome_info = "~/public/genomes/UCSC/mm10/STAR_mm10_index/chrNameLength.txt",
    log: "logs/nucleosome_pos/{group}.log"
    shell:
        '''
        java -Xmx95g -jar /home/dic/public/tools/hmmratac-1.2.10/share/hmmratac-1.2.10-1/HMMRATAC.jar -b {input.merged_bam} \
        -i {input.bam_index} -g {params.genome_info} -o nucleosome_pos/"{wildcards.group}" --bedgraph TRUE &>{log}
        '''