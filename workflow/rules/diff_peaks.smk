rule diff_peaks:
    input: 
        featureCount = "../results/{data}.featureCount.tsv",
        sample_table = "sample_table_{data}.tsv"
    output: 
        diff_peaks_report = "diff_peaks/{data}_report.html"
    threads: 1
    params:
        mem = '20G',
        jobName = "diff_peaks.{data}"
    log: "logs/diff_peaks/{data}.log"
    script:
        "../scripts/differential_peaks.Rmd"
        