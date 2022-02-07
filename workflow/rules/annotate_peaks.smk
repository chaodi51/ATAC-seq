rule annotate_peaks:
    input: 
        merged_peaks = "macs2_merged_bam/narrowPeak/{data}.merged_peaks.txt"
    output: 
        annotated_peaks = "macs2_merged_bam/narrowPeak/{data}.annotated_peaks.txt",
        anno_stats = "macs2_merged_bam/narrowPeak/{data}.annotated_peaks.stats" 
    threads: 1
    params:
        mem = '10G',
        jobName = "annotate_peaks.{data}",
        genome = "mm10"
    log: "logs/annotate_peaks/{data}.log"
    shell:
        '''
        annotatePeaks.pl {input.merged_peaks} {params.genome} -annStats {output.anno_stats}> {output.annotated_peaks} 2>> {log}    
        '''