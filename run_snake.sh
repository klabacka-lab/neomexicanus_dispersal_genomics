#!/bin/bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -J "neomex_snakemake"
#SBATCH --output=logs/snakemake.out
#SBATCH --error=logs/snakemake.err
#SBATCH --mail-user=vanwper@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#load miniforge and then activate the correct mamba environment for this project
module load miniforge3
module load apptainer
mamba activate snakemake


run the Snakefile
snakemake --jobs 2 \
  --use-conda \
  --latency-wait 60 \
  --cores 16 \
  --rerun-incomplete
