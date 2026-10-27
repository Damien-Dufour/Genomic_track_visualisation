# track visualisation - plot genomic tracks using R 

This snakemake pipeline allows to reproducibly plot genomic regions in R including 3D interactions. 

It takes as input bigwig files, bepde file in case of interaction and supports mouse and human genomes. 

# quick start 

Make sure that snakemake and R are installed/loaded 

modify __input.csv__ to give the paths of the bigwig files as well as their metadata. You can add many and select them using the config file.

modify __config.yaml__ to :
        - select which bigwig files to plot based on their modality, condition. 
        - adapt plotting
interaction file should be bedpe without header and have the columns arranged in the following order "chr1", "start1", "end1", "chr2", "start2", "end2", "count" 

Example files are provided in data/example

# run the pipeline

```bash

cd track_visualisation

snakemake --cores 1

```
