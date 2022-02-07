# merge bam files of replicates from the same condition
rule merge_bams:
    input: 
        bam = lambda wildcards: expand(
            "STAR_align_filter/{sample}.rmBList.bam", 
            sample = sample_table[sample_table["group"] == wildcards.group]["sample"].tolist())
    output: 
        merged_bam = "merged_bams/{group}.bam"
    threads: 8
    params:
        mem = '10G',
        jobName = "merge_bams.{group}"
    log: "logs/merge_bams/{group}.log"
    shell:
        '''
        # samtools merge --threads {threads} {output.merged_bam} {input.bam}  &>> {log}
        samtools index {output.merged_bam} &>> {log}
        '''