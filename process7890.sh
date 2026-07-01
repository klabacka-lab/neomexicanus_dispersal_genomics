#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G  # memory per CPU core
#SBATCH --output=logs/prcs7890.out
#SBATCH --error=logs/prcs7890.err
#SBATCH --mail-user=vanwper@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

WORKDIR="/home/vanwper/nobackup/autodelete/neomex"

module load miniforge3
module load apptainer 

mamba activate pacbioProcessing

#pbmm2 align \
	#"$WORKDIR/xeng_data/k25/7890/7890-host.fq.gz" \
	#"$WORKDIR/reference/BaumannLab/a_arizonae_AspAri2.0.mmi" \
	#"$WORKDIR/sample_bams/7890/7890-host.bam" \
	#--preset CCS \
	#--sort
#echo "host complete"

#pbmm2 align \
	#"$WORKDIR/xeng_data/k25/7890/7890-graft.fq.gz" \
	#"$WORKDIR/reference/BaumannLab/a_marmoratus_AspMarm2.0.mmi" \
	#"$WORKDIR/sample_bams/7890/7890-graft.bam" \
	#--preset CCS \
	#--sort
#echo "graft complete"

#apptainer exec \
	#--cleanenv \
	#--bind /grphome,/nobackup \
	#./docker/deepvariant.sif \
	#run_deepvariant \
		#--model_type=PACBIO \
		#--ref="$WORKDIR/reference/BaumannLab/a_arizonae_AspAri2.0.fasta" \
		#--reads="$WORKDIR/sample_bams/7890/7890-host.bam" \
		#--output_vcf="$WORKDIR/sample_vcfs/7890/7890-host.vcf.gz" \
		#--output_gvcf="$WORKDIR/sample_vcfs/7890/7890-host.g.vcf.gz" \
		#--num_shards=8
#echo "host processed"

#apptainer exec \
	#--cleanenv \
	#--bind /grphome,/nobackup \
	#./docker/deepvariant.sif \
	#run_deepvariant \
		#--model_type=PACBIO \
		#--ref="$WORKDIR/reference/BaumannLab/a_marmoratus_AspMarm2.0.fasta" \
		#--reads="$WORKDIR/sample_bams/7890/7890-graft.bam" \
		#--output_vcf="$WORKDIR/sample_vcfs/7890/7890-graft.vcf.gz" \
	       	#--output_gvcf="$WORKDIR/sample_vcfs/7890/7890-graft.g.vcf.gz" \
		#--num_shards=8
#echo "graft processed"

mamba deactivate
mamba activate bcftools

bcftools view \
	-e 'QUAL<=20' \
	"$WORKDIR/sample_vcfs/7890/7890-host.vcf.gz" \
	-Oz -o "$WORKDIR/sample_vcfs/7890/7890-host.step1.qual.vcf.gz"
bcftools index -t "$WORKDIR/sample_vcfs/7890/7890-host.step1.qual.vcf.gz"

bcftools view \
	-e 'QUAL<=20' \
	"$WORKDIR/sample_vcfs/7890/7890-graft.vcf.gz" \
	-Oz -o "$WORKDIR/sample_vcfs/7890/7890-graft.step1.qual.vcf.gz"
bcftools index -t "$WORKDIR/sample_vcfs/7890/7890-graft.step1.qual.vcf.gz"

bcftools view \
	-i 'FMT/GQ>=20' \
	"$WORKDIR/sample_vcfs/7890/7890-host.step1.qual.vcf.gz" \
	-Oz -o "$WORKDIR/sample_vcfs/7890/7890-host.step2.gq.vcf.gz"
bcftools index -t "$WORKDIR/sample_vcfs/7890/7890-host.step2.gq.vcf.gz"

bcftools view \
	-i 'FMT/GQ>=20' \
	"$WORKDIR/sample_vcfs/7890/7890-graft.step1.qual.vcf.gz" \
	-Oz -o "$WORKDIR/sample_vcfs/7890/7890-graft.step2.gq.vcf.gz"
bcftools index -t "$WORKDIR/sample_vcfs/7890/7890-graft.step2.gq.vcf.gz"

bcftools view \
	-i 'FMT/DP>=10' \
	"$WORKDIR/sample_vcfs/7890/7890-graft.step2.gq.vcf.gz" \
	-Oz -o "$WORKDIR/sample_vcfs/7890/7890-graft.step3.dp.vcf.gz"
bcftools index -t "$WORKDIR/sample_vcfs/7890/7890-graft.step3.dp.vcf.gz"

bcftools view \
	-i 'FMT/DP>=10' \
	"$WORKDIR/sample_vcfs/7890/7890-host.step2.gq.vcf.gz" \
	-Oz -o "$WORKDIR/sample_vcfs/7890/7890-host.step3.dp.vcf.gz"
bcftools index -t "$WORKDIR/sample_vcfs/7890/7890-host.step3.dp.vcf.gz"

bcftools view \
	-i 'GT="1/1"' \
	"$WORKDIR/sample_vcfs/7890/7890-host.step3.dp.vcf.gz" \
	-Oz -o "$WORKDIR/sample_vcfs/7890/7890-host.step4.gt.vcf.gz"
bcftools index -t "$WORKDIR/sample_vcfs/7890/7890-host.step4.gt.vcf.gz"

bcftools view \
	-i 'GT="1/1"' \
	"$WORKDIR/sample_vcfs/7890/7890-graft.step3.dp.vcf.gz" \
	-Oz -o "$WORKDIR/sample_vcfs/7890/7890-graft.step4.gt.vcf.gz"
bcftools index -t "$WORKDIR/sample_vcfs/7890/7890-graft.step4.gt.vcf.gz"
