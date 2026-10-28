# track visualisation - plot genomic tracks using R 

A [Snakemake](https://snakemake.readthedocs.io/) tool that plots genomic tracks from bigWig files in a genomic window defined by the user, optionally overlaying 3D chromatin interactions from a bedpe file. Supports mouse and human genomes.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/) ≥ 7
- R ≥ 4.0 with bioconductor 

```bash
# quick check
snakemake --version
Rscript -e 'packageVersion("BiocManager")' 
```

## Directory structure

```
track_visualisation/
├── Snakefile                # workflow definition
├── config
├   └── config.yaml              # default configuration
├── data/
├   ├── input.csv            # sample sheet (all candidate bigWigs + metadata)
├   └── example/             # small example dataset for a test run
├── results/                 # output folder
└── track_visualisation.R    # plotting script
```

## Inputs

### Sample sheet — `data/input.csv`

The pipeline selects the bigWig files to plot from this sheet, filtered by the conditions set in the config file. Expected columns:

| column      | description                              |
|-------------|------------------------------------------|
| `path`      | path to the bigWig files                  |
| `condition` | treatment, e.g. `vehicle` or `drug`       |
| `day`       | time point, e.g. `0`                     |
| `modality`  | assay, e.g. `ATAC`                       |

Example:

```csv
path,condition,day,modality
data/example/sample1.bw,vehicle,0,ATAC
data/example/sample2.bw,drug,0,ATAC
```

### Interaction file (optional) — bedpe, no header

Columns must be in this exact order:
`chr1, start1, end1, chr2, start2, end2, count`

```bash
chr3	12316965	12317413	chr3	12431759	12449688	217
chr3	12270023	12310329	chr3	12429070	12446143	552
chr3	12319127	12386784	chr3	12293468	12333145	292
```

Example files are provided in `data/example/`.

## Configuration

All parameters live in `config.yaml`:

| key                | type    | description                                            | example |
|--------------------|---------|--------------------------------------------------------|---------|
| `title`            | string  | analysis name, used as output subfolder                | `Gapdh_locus` |
| `bw_file`          | path    | sample sheet                                           | `data/input.csv` |
| `interaction_file` | path    | bedpe file with 3D interactions                         | `data/example/x.bedpe` |
| `results_dir`      | path    | output directory                                        | `results/` |
| `coordinates`      | string  | genomic window, format `chr:start-end`                  | `"chr3:12,252,751-12,510,000"` |
| `max`              | number  | y-axis maximum per modality                             | `7` |
| `day`              | number  | which day(s) to plot                                    | `"day0,day3"` |
| `modality`         | string  | which modality to plot                                  | `ATAC` |
| `ttt`              | list    | which treatments to include                             | `["vehicle, drug"]` |
| `interaction`      | boolean | overlay the interaction file (`TRUE`/`FALSE`)           | `TRUE` |

Keeping one config per analysis is preferred (e.g. `config/Gapdh_locus_2026-10-27.yaml`).
Rerunning with an old config reproduces the exact figure — a built-in log of past plots.

## Usage

```bash
cd track_visualisation
snakemake --cores 1
```

Useful options:

```bash
snakemake -n                 # dry run: shows what would be done
snakemake -p                 # print the executed shell commands
snakemake --configfile config/my_analysis.yaml   # use a different config
```

## Outputs

Written to `<results_dir>/<title>/`:

```
results/
└── Gapdh_locus/
    ├── figures/Gapdh_locus_track.pdf    # final track figure
    ├── logs/test_rule.log               # job log
    └── data/Gapdh_locus_input_data.csv  # (intermediate) parsed sample sheet
```
