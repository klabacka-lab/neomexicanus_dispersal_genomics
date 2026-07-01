#!/bin/bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --cpus-per-task=8
#SBATCH --mem=84G
#SBATCH -J "neomex_snakemake"
#SBATCH --output=logs/snakemake.out
#SBATCH --error=logs/snakemake.err
#SBATCH --mail-user=vanwper@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --requeue

#load miniforge and then activate the correct mamba environment for this project
module load miniforge3
module load apptainer
mamba activate snakemake

snakemake --unlock

#run the Snakefile
snakemake --jobs 1 \
  --use-conda \
  --latency-wait 60 \
  --cores 8 \
  --rerun-incomplete
