rule featureCount:
    input:
        bamfiles = lambda wildcards: expand(
            "STAR_align_filter/{sample}.rmBList.bam", 
            sample = sample_table[sample_table["data"] == wildcards.data]["sample"].tolist()),
        merged_peaks_saf = "macs2_merged_bam/narrowPeak/{data}.merged_peaks.saf"
    output:
        featureCount = "../results/{data}.featureCount.tsv"
    log:
        "logs/featureCount/{data}.log"   
    threads: 8     
    params: 
      mem = '10G',
      jobName = "{data}.featureCount",
      strand = 0,
      fracOverlap = 0,
      fracOverlapFeature = 0
    script:
        "../scripts/featureCount_segments.R"

## return a list of library strand, as 'none', 'yes' or 'reverse' denoted as in HTSeq -s option
# def get_strandness(sample_table):
#     if "strand" in sample_table.columns:
#         return sample_table["strand"].tolist()

# # def get_deseq2_threads(wildcards=None):
# #     # https://twitter.com/mikelove/status/918770188568363008
# #     few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
# #     return 1 if len(sample_list) < 100 or few_coeffs else 6   

# rule count_matrix:
#     input:
#         expand("STAR_align/{sample}.ReadsPerGene.out.tab", sample=SAMPLES)
#     output:
#         "../results/tables/all_readCount.tsv"
#     params:
#         sample_list = sample_table["sample"].tolist(),
#         strand_list = get_strandness(sample_table),
#         mem = '2G',
#         jobName = "count_matrix" 
#     script:
#         "../scripts/count-matrix.py"

