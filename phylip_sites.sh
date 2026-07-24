#!/bin/bash

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16G   # memory per CPU core
#SBATCH --mail-user=vanwper@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/sites.out
#SBATCH --error=logs/sites.err

workingdir="/home/vanwper/nobackup/autodelete/neomex/"

module load miniforge3
mamba activate python

python vcf2phylip.py \
	-i "$workingdir/analysis/graft-final-SNPs.vcf.gz" \
	--output-folder "$workingdir/misc" \
	--output-prefix "graft-SNPs" \
	--fasta \
	-w

python vcf2phylip.py \
        -i "$workingdir/analysis/host-final-SNPs.vcf.gz" \
        --output-folder "$workingdir/misc" \
        --output-prefix "host-SNPs" \
        --fasta \
        -w
