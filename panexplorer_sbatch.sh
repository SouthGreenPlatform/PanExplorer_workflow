#!/bin/bash

#SBATCH --job-name=panexplorer
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=supermem

module load singularity/4.0.1
export PANEX_PATH=$PWD

singularity exec $PANEX_PATH/singularity/panexplorer.sif snakemake --cores 1 -s Snakemake_files/Snakefile_pggb_heatmap_upset
