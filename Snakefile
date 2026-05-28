#define locations for the parent reference fastas
P1_REF = "reference/BaumannLab/a_arizonae_AspAri2.0.fasta.gz"
P2_REF = "reference/BaumannLab/a_marmoratus_AspMarm2.0.fasta.gz"

#sample numbers
IDS = ["7890","7891","7892","7893","7894","7895","7896","7897"]

#xengsort classify output types
XENG_CLASS = ["host", "graft", "both", "neither", "ambiguous"]

#varying kmer lengths to determine optimal kmer length
KMER = [20, 23, 25, 28, 30]

rule all:
  input:
    expand("sample_bams/{sample}/{sample}-host.bam", sample=IDS),
    expand("sample_bams/{sample}/{sample}-graft.bam", sample=IDS)

rule xengsort_index:
  input:
    host=P1_REF,
    graft=P2_REF
  output:
      "xeng_data/k{k}/neomex.hash",
      "xeng_data/k{k}/neomex.info"
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
      hash_file="xeng_data/k{k}/neomex.hash",
      fastq="sample_fastqs/{ids}/{ids}.fq"
  output:
      "xeng_data/k{k}/{ids}/{ids}-{xeng}.fq.gz"
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
    "reference/BaumannLab/{name}.fasta"
  output:
    "reference/BaumannLab/{name}.mmi"
  conda:
    "/home/vanwper/.conda/envs/pacbioProcessing"
  shell:
    """
    pbmm2 index {input} {output}
    """

rule pbmm2_align_host:
  input:
    reads="xeng_data/k25/{sample}/{sample}-host.fq.gz",
    ref="reference/BaumannLab/a_arizonae_AspAri2.0.mmi"
  output:
    "sample_bams/{sample}/{sample}-host.bam"
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
    reads="xeng_data/k25/{sample}/{sample}-graft.fq.gz",
    ref="reference/BaumannLab/a_marmoratus_AspMarm2.0.mmi"
  output:
    "sample_bams/{sample}/{sample}-graft.bam"
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
