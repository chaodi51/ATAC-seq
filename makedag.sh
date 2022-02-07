
snakemake --rulegraph | dot -Tsvg > DAG.svg
rsvg-convert DAG.svg > DAG.png
