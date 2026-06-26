#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G  # memory per CPU core
#SBATCH --output=logs/prcs7890.out
#SBATCH --error=logs/prcs7890.err
#SBATCH --mail-user=vanwper@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

WORKDIR="/home/vanwper/nobackup/autodelete/neomex"

module load miniforge3

mamba activate pacbioProcessing

pbmm2 align \
	"$WORKDIR/xeng_data/k25/7890/7890-host.fq.gz" \
	"$WORKDIR/reference/BaumannLab/a_arizonae_AspAri2.0.mmi" \
	"$WORKDIR/sample_bams/7890/7890-host.bam" \
	--preset CCS \
	--sort
echo "host complete"

pbmm2 align \
	"$WORKDIR/xeng_data/k25/7890/7890-graft.fq.gz" \
	"$WORKDIR/reference/BaumannLab/a_marmoratus_AspMarm2.0.mmi" \
	"$WORKDIR/sample_bams/7890/7890-graft.bam" \
	--preset CCS \
	--sort
echo "graft complete"
