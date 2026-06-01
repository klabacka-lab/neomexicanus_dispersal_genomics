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

#convert unfiltered bam file into a fastq
#mamba activate samtools_env
#samtools fastq $bam | bgzip -c > $workingdir/7888/7888.fastq.gz

#use xengsort to classify the reads of the sample by which parent it comes from host is arizonae and graft is marmoratus
#mamba activate xengsort
#xengsort classify \
       	#--index /home/vanwper/nobackup/autodelete/neomex/xeng_data/k25/neomex.hash \
       	#--fastq /home/vanwper/nobackup/autodelete/neomex/7888/7888.fastq.gz \
	#--prefix /home/vanwper/nobackup/autodelete/neomex/7888/7888 \
	#--mode count \
	#--compression gz

#Index for the reference files is already created so just need to align
#the graft and host files to their respective index
mamba activate pacbioProcessing
pbmm2 align \
	$workingdir/reference/BaumannLab/a_arizonae_AspAri2.0.mmi \
	$workingdir/7888/7888-host.fq.gz \
	$workingdir/7888/7888-host.bam \
	--preset CCS \
	--sort
echo "host processing complete"

pbmm2 align \
	$workingdir/reference/BaumannLab/a_marmoratus_AspMarm2.0.mmi \
	$workingdir/7888/7888-graft.fq.gz \
	$workingdir/7888/7888-graft.bam \
	--preset CCS \
	--sort
echo "processing complete"
