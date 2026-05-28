#!/bin/bash

#SBATCH --time=8:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16384M   # memory per CPU core
#SBATCH --output=logs/prcs7888.out
#SBATCH --error=logs/prcs7888.err
#SBATCH --mail-user=vanwper@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

workingdir="/home/vanwper/nobackup/autodelete/neomex"
bam="/home/vanwper/groups/grp_reptile_genomics/nobackup/archive/raw_neomex_genomic_data/7888/m84100_250712_050951_s4.hifi_reads.bc1030.bam"

module load miniforge3
mamba activate samtools_env

samtools fastq $bam | bgzip -c > $workingdir/7888/7888.fastq.gz
