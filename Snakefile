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
    f"{WORKDIR}/analysis/merged-graft.bcf"

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

#filters vcf files by QUAL score
rule variant_filtration_step1:
  input:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.vcf.gz"
  output:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.step1.qual.vcf.gz"
  params:
    qual=20
  conda:
    "/home/vanwper/.conda/envs/bcftools"
  shell:
    """
    bcftools view \
      -e 'QUAL<={params.qual}' \
      {input.vcf} \
      -Oz -o {output.vcf}

    bcftools index -t {output.vcf}
    """

#filters the existing QUAL filtered vcf files by GQ score
rule variant_filtration_step2:
  input:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.step1.qual.vcf.gz"
  output:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.step2.gq.vcf.gz"
  params:
    gq=20
  conda:
    "/home/vanwper/.conda/envs/bcftools"
  shell:
    """
    bcftools view \
      -i 'FMT/GQ>={params.gq}' \
      {input.vcf} \
      -Oz -o {output.vcf}

    bcftools index {output.vcf}
    """

#filters GQ filtered files by read depth
rule variant_filtration_step3:
  input:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.step2.gq.vcf.gz"
  output:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.step3.dp.vcf.gz"
  params:
    dp=10
  conda:
    "/home/vanwper/.conda/envs/bcftools"
  shell:
    """
    bcftools view \
      -i 'FMT/DP>={params.dp}' \
      {input.vcf} \
      -Oz -o {output.vcf}

    bcftools index {output.vcf}
    """

#filters GT data
rule variant_filtration_step4:
  input:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.step3.dp.vcf.gz"
  output:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.step4.gt.vcf.gz"
  conda:
    "/home/vanwper/.conda/envs/bcftools"
  shell:
    """
    bcftools view \
      -i 'GT="1/1"' \
      {input.vcf} \
      -Oz -o {output.vcf}
  
    bcftools index {output.vcf}
    """

#creates a bed file to determine the sites of interest
rule union_sites_host:
  input:
    expand(
      f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-host.step4.gt.vcf.gz",
      sample=IDSP
    )
  output:
    f"{WORKDIR}/analysis/union_sites_host.bed"
  conda:
    "/home/vanwper/.conda/envs/bcftools"
  shell:
    """
    bcftools query -f '%CHROM\\t%POS0\\t%POS\\n' {input} \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > {output}
    """

rule union_sites_graft:
  input:
    expand(
      f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-graft.step4.gt.vcf.gz",
      sample=IDSP
    )
  output:
    f"{WORKDIR}/analysis/union_sites_graft.bed"
  conda:
    "/home/vanwper/.conda/envs/bcftools"
  shell:
    """
    bcftools query -f '%CHROM\\t%POS0\\t%POS\\n' {input} \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > {output}
    """

#create gvcf to determine if data is missing
rule site_specific_variant_call_host:
  input:
    bam=f"{WORKDIR}/sample_bams/{{sample}}/{{sample}}-host.bam",
    ref=f"{WORKDIR}/reference/BaumannLab/a_arizonae_AspAri2.0.fasta"
  threads: 64
  output:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-host.dupe.vcf.gz",
    gvcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-host.g.vcf.gz"
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
        --output_gvcf={output.gvcf} \
        --num_shards={threads}
    """

rule site_specific_variant_call_graft:
  input:
    bam=f"{WORKDIR}/sample_bams/{{sample}}/{{sample}}-graft.bam",
    ref=f"{WORKDIR}/reference/BaumannLab/a_marmoratus_AspMarm2.0.fasta"
  threads: 64
  output:
    vcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-graft.dupe.vcf.gz",
    gvcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-graft.g.vcf.gz"
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
        --output_gvcf={output.gvcf} \
        --num_shards={threads}
    """

rule name_samples:
  input:
    gvcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.g.vcf.gz"
  output:
    gvcf=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.fixed.g.vcf.gz",
    tbi=f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-{{xeng}}.fixed.g.vcf.gz.tbi"
  params:
    sample=lambda wildcards: wildcards.sample
  conda:
    "/home/vanwper/.conda/envs/bcftools"
  shell:
    """
    sample_file=$(mktemp)

    echo "{params.sample}" > "$sample_file"

    bcftools reheader \
      -s "$sample_file" \
      {input.gvcf} | \
    bcftools view -Oz -o {output.gvcf}

    bcftools index -f -t {output.gvcf}

    rm -f "$sample_file"
    """

rule merge_host_gvcf:
  input:
    gvcfs=expand(f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-host.fixed.g.vcf.gz", sample=IDSP)
  output:
    bcf=f"{WORKDIR}/analysis/merged-host.bcf"
  conda:
    "/home/vanwper/.conda/envs/pacbioProcessing"
  shell:
    """
    rm -rf GLnexus.DB

    apptainer exec \
      --cleanenv \
      --bind /grphome,/nobackup \
      ./docker/glnexus.sif \
      glnexus_cli \
        --config DeepVariant \
        --threads 4 \
        --mem-gbytes 48 \
        {input.gvcfs} \
        > {output.bcf}
    """

rule merge_graft_gvcf:
  input:
    gvcfs=expand(f"{WORKDIR}/sample_vcfs/{{sample}}/{{sample}}-graft.fixed.g.vcf.gz", sample=IDSP)
  output:
    bcf=f"{WORKDIR}/analysis/merged-graft.bcf"
  conda:
    "/home/vanwper/.conda/envs/pacbioProcessing"
  shell:
    """
    rm -rf GLnexus.DB

    apptainer exec \
      --cleanenv \
      --bind /grphome,/nobackup \
      ./docker/glnexus.sif \
      glnexus_cli \
        --config DeepVariant \
        --threads 4 \
        --mem-gbytes 48 \
        {input.gvcfs} \
        > {output.bcf}
    """
