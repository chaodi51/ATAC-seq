# merge bam files of replicates from the same condition
rule footprinting:
    input: 
        merged_bam = "merged_bams/{group}.bam",
        peak = "macs2_merged_bam/narrowPeak/{group}_peaks.narrowPeak"
    output: 
        # motif = "footprinting/{group}.bed",
        wig = "footprinting/{group}.wig"
    threads: 1
    params:
        mem = '40G',
        jobName = "footprinting.{group}",
        group_name = "{group}",
    log: "logs/footprinting/{group}.log"
    shell:
        '''
        # rgt-hint footprinting --atac-seq --paired-end --organism=mm10 \
        --output-location=footprinting --output-prefix={wildcards.group} {input.merged_bam} {input.peak} &> {log}
        
        rgt-hint tracks --bc --bigWig --organism=mm10 --output-location=footprinting --output-prefix={wildcards.group} {input.merged_bam} {input.peak} &>> {log}

        rgt-motifanalysis matching --organism=mm10 --output-location=footprinting --input-files footprinting/{params.group_name}.bed &>> {log}
        '''

## TBD
# rule diff_TF:
#     input:
#         matched_motif_dsbAB = "footprinting/ctrl_mpbs.bed, footprinting/dsb_mpbs.bed",
#         matched_motif_DSB = "footprinting/no_DSB_mpbs.bed, footprinting/with_DSB_mpbs.bed",
#         bam_dsbAB = "merged_bams/ctrl.bam,merged_bams/dsb.bam",
#         bam_DSB = "merged_bams/no_DSB.bam,merged_bams/with_DSB.bam",          
#     output:
#         
#     threads:8
#     params:
#         mem = '40G',
#         jobName = "diff_TF"
#     log: "logs/diff_TF.log"
#     shell:     
#         '''
#         rgt-hint differential --organism=mm10 --bc --nc {threads} --mpbs-files={input.matched_motif_dsbAB} \
#          --reads-files={input.bam_dsbAB} --conditions=ctrl,dsb --output-location=footprinting/diff_dsbAB &>> {log}
# 
#       rgt-hint differential --organism=mm10 --bc --nc {threads} --mpbs-files={input.matched_motif_DSB} \
#          --reads-files={input.bam_DSB} --conditions=no_DSB,with_DSB --output-location=footprinting/diff_DSB &>> {log} 
        # '''
