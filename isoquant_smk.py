SAMPLES = ["TM-1_ovule_0dpa-1_L2_368X68", "TM-1_ovule_0dpa-2_L2_369X69", "TM-1_ovule_0dpa-3_L1_370X70", "TM-1_fiber_10dpa-1_L2_380X80", "TM-1_fiber_10dpa-2_L2_381X81", "TM-1_fiber_10dpa-3_L2_383X83", "TM-1_fiber_15dpa-1_L2_384X84", "TM-1_fiber_15dpa-2_L2_386X86", "TM-1_fiber_15dpa-3_L3_387X87", "TM-1_fiber_20dpa-1_L2_354X54", "TM-1_fiber_20dpa-2_L2_389X89", "TM-1_fiber_20dpa-3_L2_355X55", "TM-1_fiber_25dpa-1_L2_391X91", "TM-1_fiber_25dpa-2_L2_392X92", "TM-1_fiber_25dpa-3_L2_362X62"]

rule all:
  input:
    ["clean_data/{}.R{}_fp.fastq.gz".format(sample, i) for sample in SAMPLES for i in [1, 2]],
    expand("quantifying/{sample}", sample = SAMPLES),
    matrix = "Mergequantifying/Salmon"

#使用fastp处理数据
rule trim_fastq:
  input:
    r1="raw_data/{sample}.R1.fastq.gz",
    r2="raw_data/{sample}.R2.fastq.gz"
  output:
    r1="clean_data/{sample}.R1_fp.fastq.gz",
    r2="clean_data/{sample}.R2_fp.fastq.gz",
    report="clean_data/{sample}_fp.html",
    json="clean_data/{sample}_fp.json"
  shell:
    """
      fastp -i {input.r1} -I {input.r2} \
      -o {output.r1} -O {output.r2} \
      -h {output.report} -j {output.json}
    """

#Salmon无比对定量
#构建索引
rule salmon_index:
  input:
    transcriptome="Ghirsutum_transcript.fa"
  output:
    "transcripts"
  shell:
    """
      salmon index -t {input.transcriptome} -i transcripts --keepDuplicates
    """
#定量
rule quant_salmon:
  input:
    r1="clean_data/{sample}.R1_fp.fastq.gz",
    r2="clean_data/{sample}.R2_fp.fastq.gz",
    idx="transcripts"
  output:
    quanti= "quantifying/{sample}"
  shell:
    """
      salmon quant -i {input.idx} -l A -1 {input.r1} \
      -2 {input.r2} -o {output.quanti}
    """
     
#合并表达矩阵
rule abundance_estimates_to_matrix:
  input:
    counts = expand("quantifying/{sample}", sample=SAMPLES)
  output:
    matrix = "Mergequantifying/Salmon"
  shell:
    """
      for sample in {input.counts}; do
        echo $sample/quant.sf >> Mergequantifying/quant_files.fofn
      done
      awk '{{print $10"\t"$12}}' Ghincfib.gtf | sed 's/"//g' | sed 's/;//g' | sort | uniq | sed '1d' > Mergequantifying/gene_trans_map
      perl /public/home/software/opt/bio/software/Trinity/2.11.0/util/abundance_estimates_to_matrix.pl --est_method salmon \
      --quant_files Mergequantifying/quant_files.fofn --name_sample_by_basedir --gene_trans_map Mergequantifying/gene_trans_map \
      --out_prefix {output.matrix}
      sed '1s/^/isoform_id/' Mergequantifying/Salmon.isoform.TMM.EXPR.matrix > Mergequantifying/Ghir.isoform.TMM.EXPR.matrix
    """



