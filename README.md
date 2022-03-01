#Snakemake workflow: ATAC-seq

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/ATAC-seq.svg?branch=master)](https://travis-ci.org/snakemake-workflows/ATAC-seq)

This is using the standard Snakemake workflow template. Replace this text with a comprehensive description covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules`, and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

This is the first version of ATAC-seq, using Bassing lab datasets dsb vs ctrl and no_DSB vs with_DSB
## Authors

* chaodi (dic@chop.edu)

## Usage
Running on new respublica by:
snakemake --latency-wait 10 -j 10 -p -c "sbatch --job-name={params.jobName} --mem={params.mem} -c {threads} --time=360 -e sbatch/{params.jobName}.e -o sbatch/{params.jobName}.o"

## Workflow
![alt text](https://github.com/chaodi51/ATAC-seq/blob/master/workflow/DAG.png?raw=true)
