#!/bin/bash

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32768M   # memory per CPU core
#SBATCH --mail-user=vanwper@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/index_ref.out
#SBATCH --error=logs/index_ref.err

workingdir="/home/vanwper/nobackup/autodelete/neomex/reference/BaumannLab"

module load miniforge3

#use samtools to index the reference fastas
mamba activate samtools_env

samtools faidx $workingdir/a_arizonae_AspAri2.0.fasta
echo "arizonae complete"

samtools faidx $workingdir/a_marmoratus_AspMarm2.0.fasta
echo "marmoratus complete"
