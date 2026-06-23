#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8192M   # memory per CPU core
WORKDIR="/home/vanwper/nobackup/autodelete/neomex"

module load miniforge3

mamba activate pacbioProcessing

for f in "$WORKDIR"/sample_vcfs/*/*-graft.fixed.g.vcf.gz; do
	echo "Checking $f"
	gatk ValidateVariants \
	       	-V "$f" \
	       	-R "$WORKDIR"/reference/BaumannLab/a_marmoratus_AspMarm2.0.fasta	
done
