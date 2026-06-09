#define locations for the parent reference fastas
P1_REF = "{WORKDIR}/reference/BaumannLab/a_arizonae_AspAri2.0.fasta.gz"
P2_REF = "{WORKDIR}/reference/BaumannLab/a_marmoratus_AspMarm2.0.fasta.gz"

#define working directory for files
WORKDIR = "/home/vanwper/nobackup/autodelete/neomex"

#sample numbers
IDS = ["7891","7892","7893","7894","7895","7896","7897"]
#sample numbers including the updated 7888
IDSP = ["7888", "7891","7892","7893","7894","7895","7896","7897"]

#xengsort classify output types
XENG_CLASS = ["host", "graft", "both", "neither", "ambiguous"]

#varying kmer lengths to determine optimal kmer length
KMER = [20, 23, 25, 28, 30]

rule all:
  input:
    expand(f"{WORKDIR}/sample_vcfs/{{id}}/{{id}}-host.vcf.gz", id=IDSP),
    expand(f"{WORKDIR}/sample_vcfs/{{id}}/{{id}}-graft.vcf.gz", id=IDSP)

rule xengsort_index:
  input:
    host=P1_REF,
    graft=P2_REF
  output:
      "{WORKDIR}/xeng_data/k{k}/neomex.hash",
      "{WORKDIR}xeng_data/k{k}/neomex.info"
  conda:
      "/home/vanwper/.conda/envs/xengsort"
  shell:
      "xengsort index --index xeng_data/k{wildcards.k}/neomex "
      " -H {input.host} "
      " -G {input.graft} "
      "-k {wildcards.k} "
      "-n 2_500_000_000"


rule xengsort_classify:
  input:
      hash_file="{WORKDIR}/xeng_data/k{k}/neomex.hash",
      fastq="{WORKDIR}/sample_fastqs/{ids}/{ids}.fq"
  output:
      "{WORKDIR}/xeng_data/k{k}/{ids}/{ids}-{xeng}.fq.gz"
  params:
      mode="count"
  conda:
      "/home/vanwper/.conda/envs/xengsort"
  shell:
      """
      xengsort classify \
          --index {input.hash_file} \
          --fastq {input.fastq} \
          --prefix xeng_data/k{wildcards.k}/{wildcards.ids}/{wildcards.ids} \
          --classification {params.mode} \
          --compression gz
      """

rule pbmm2_index:
  input:
    "{WORKDIR}/reference/BaumannLab/{name}.fasta"
  output:
    "{WORKDIR}/reference/BaumannLab/{name}.mmi"
  conda:
    "/home/vanwper/.conda/envs/pacbioProcessing"
  shell:
    """
    pbmm2 index {input} {output}
    """

rule pbmm2_align_host:
  input:
    reads="{WORKDIR}/xeng_data/k25/{sample}/{sample}-host.fq.gz",
    ref="{WORKDIR}/reference/BaumannLab/a_arizonae_AspAri2.0.mmi"
  output:
    "{WORKDIR}/sample_bams/{sample}/{sample}-host.bam"
  conda:
    "/home/vanwper/.conda/envs/pacbioProcessing"
  shell:
    """
    pbmm2 align \
      {input.ref} \
      {input.reads} \
      {output} \
      --preset CCS \
      --sort
    """

rule pbmm2_align_graft:
  input:
    reads="{WORKDIR}/xeng_data/k25/{sample}/{sample}-graft.fq.gz",
    ref="{WORKDIR}/reference/BaumannLab/a_marmoratus_AspMarm2.0.mmi"
  output:
    "{WORKDIR}/sample_bams/{sample}/{sample}-graft.bam"
  conda:
    "/home/vanwper/.conda/envs/pacbioProcessing"
  shell:
    """
    pbmm2 align \
      {input.ref} \
      {input.reads} \
      {output} \
      --preset CCS \
      --sort
    """

rule variant_call_host:
  input:
    bam=f"{WORKDIR}/sample_bams/{{sample}}/{{sample}}-host.bam",
    ref=f"{WORKDIR}/reference/BaumannLab/a_arizonae_AspAri2.0.fasta"
  threads: 8
  output:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-host.vcf.gz"
  shell:
    """
    apptainer exec \
      --cleanenv \
      --bind /grphome,/nobackup \
      ./docker/deepvariant.sif \
      run_deepvariant \
        --model_type=PACBIO \
        --ref={input.ref} \
        --reads={input.bam} \
        --output_vcf={output.vcf} \
        --num_shards={threads}
    """

rule variant_call_graft:
  input:
    bam=f"{WORKDIR}/sample_bams/{{sample}}/{{sample}}-graft.bam",
    ref=f"{WORKDIR}/reference/BaumannLab/a_marmoratus_AspMarm2.0.fasta"
  threads: 8
  output:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-graft.vcf.gz"
  shell:
    """
    apptainer exec \
      --cleanenv \
      --bind /grphome,/nobackup \
      ./docker/deepvariant.sif \
      run_deepvariant \
        --model_type=PACBIO \
        --ref={input.ref} \
        --reads={input.bam} \
        --output_vcf={output.vcf} \
        --num_shards={threads}
    """
