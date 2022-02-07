rule merge_peaks:
    input: 
        allpeaks = lambda wildcards: expand(
            "macs2_merged_bam/narrowPeak/{group}_peaks.narrowPeak", 
            group = sample_table[sample_table["data"] == wildcards.data]["group"].tolist())
    output: 
        merged_peaks = "macs2_merged_bam/narrowPeak/{data}.merged_peaks.txt",
        merged_peaks_saf = "macs2_merged_bam/narrowPeak/{data}.merged_peaks.saf"
    threads: 1
    params:
        mem = '5G',
        jobName = "merge_peaks.{data}",
    log: "logs/merge_peaks/{data}.log"
    shell:
        '''
        mergePeaks {input.allpeaks} > {output.merged_peaks} 2>> {log}   
        sed 1d {output.merged_peaks} | awk 'BEGIN{{FS=OFS="\t"}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' > {output.merged_peaks_saf} 2>> {log}  
        '''