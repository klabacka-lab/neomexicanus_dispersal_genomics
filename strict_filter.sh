#!/bin/bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -J "neomex_snakemake"
#SBATCH --output=logs/strict_filters.out
#SBATCH --error=logs/strict_filters.err
#SBATCH --mail-user=vanwper@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

WORKDIR="/home/vanwper/nobackup/autodelete/neomex"

module load miniforge3

#FILTER BY MISSING
mamba activate bcftools
bcftools view \
    -i 'F_MISSING=0' \
    "$WORKDIR/analysis/merged-filtered-graft.bcf" \
    -Ob \
    -o "$WORKDIR/strict_filters/graft-final.bcf"

bcftools index -f $WORKDIR/strict_filters/graft-final.bcf

bcftools view \
    -i 'F_MISSING=0' \
    "$WORKDIR/analysis/merged-filtered-host.bcf" \
    -Ob \
    -o "$WORKDIR/strict_filters/host-final.bcf"

bcftools index -f "$WORKDIR/strict_filters/host-final.bcf"
echo "Missingness filtration complete"


#FILTER TO SNPs
bcftools view \
    -m2 -M2 \
    -v snps \
    -Ob \
    -o "$WORKDIR/strict_filters/graft-final-SNPs.bcf" \
    "$WORKDIR/strict_filters/graft-final.bcf"

bcftools index -f "$WORKDIR/strict_filters/graft-final-SNPs.bcf"

bcftools view \
    -m2 -M2 \
    -v snps \
    -Ob \
    -o "$WORKDIR/strict_filters/host-final-SNPs.bcf" \
    "$WORKDIR/strict_filters/host-final.bcf"

bcftools index -f "$WORKDIR/strict_filters/host-final-SNPs.bcf"
echo "SNP filter complete"


#CONVERT TO VCF
bcftools view -Ov -o "$WORKDIR/strict_filters/graft-final-SNPs.vcf" "$WORKDIR/strict_filters/graft-final-SNPs.bcf"
bcftools index -t "$WORKDIR/strict_filters/graft-final-SNPs.vcf"

bcftools view -Ov -o "$WORKDIR/strict_filters/host-final-SNPs.vcf" "$WORKDIR/strict_filters/host-final-SNPs.bcf"
bcftools index -t "$WORKDIR/strict_filters/host-final-SNPs.vcf"
echo "BCFs converted to VCFs"


#CONVERT TO FASTA
mamba activate python
python vcf2phylip.py \
    -i "$WORKDIR/strict_filters/graft-final-SNPs.vcf" \
    --output-folder "$WORKDIR/strict_filters/phylogeny" \
    --output-prefix "graft-SNPs" \
    --fasta

python vcf2phylip.py \
    -i "$WORKDIR/strict_filters/host-final-SNPs.vcf" \
    --output-folder "$WORKDIR/strict_filters/phylogeny" \
    --output-prefix "host-SNPs" \
    --fasta
echo "FASTAs generated"


#GENERATE PHYLOGENY
mamba activate pacbioProcessing
iqtree3 \
    -s "$WORKDIR/strict_filters/phylogeny/graft-SNPs.min4.fasta" \
    -m MFP \
    -bb 1000 \
    -bnni \
    -T AUTO \
    --redo \
    -pre "$WORKDIR/strict_filters/phylogeny/graft-snps_tree"

iqtree3 \
    -s "$WORKDIR/strict_filters/phylogeny/host-SNPs.min4.fasta" \
    -m MFP \
    -bb 1000 \
    -bnni \
    -T AUTO \
    --redo \
    -pre "$WORKDIR/strict_filters/phylogeny/host-snps_tree"
echo "Phylogenies complete"