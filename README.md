
# 13CFluxRatioAnalyisWithFRAPPPE

All the scripts to reproduce all the analysis and figures in the paper.

Data is stored using git-lfs in this repo.

## Setup

The following tools are needed:
- pixi
- snakemake >= 8.0
- git

## Running the analysis

After placing the data in the `Data` directory, run the pipeline with the following command:
```
snakemake --workflow-profile Profiles/Slurm
```
