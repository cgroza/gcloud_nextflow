params.reads          = "reads.csv"
params.reference      = "reference.fa"
params.vcf            = "variants.vcf"
params.out            = "out"

process preprocess {
  cpus 10
  memory "60G"
  time "3h"

  input:
  path(vcf)
  path(fasta)

  output:
  path("index")

  script:
  panvcf = vcf.getBaseName()
  """
  mkdir index
  gunzip ${vcf}
  PanGenie-index -v ${panvcf} -r ${fasta} -t 10 -o index/processed
  """
}

process pangenie {
  cpus 10
  memory "60G"
  time "6h"
  publishDir "${params.out}/genotypes", mode: 'copy'

  input:
  tuple val(sample_name), path(sample_bam), path(panref), path("index/")

  output:
  path("${sample_name}_genotyping.vcf.gz")

  script:
  """
  PanGenie -t 10 -j 10 -s ${sample_name} -i <(samtools fastq --reference ${cram_ref} ${sample_bam}) -f index/index/processed -o ${sample_name}
  bgzip ${sample_name}_genotyping.vcf
  """
}

workflow {
  Channel.fromPath(params.cram_reference).set{ref_ch}
  Channel.fromPath(params.pangenome_reference).set{pan_ref_ch}
  Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path, checkIfExists:false)]}.set{reads_ch}

  Channel.fromPath(params.vcf).set{vcf_ch}

  preprocess(vcf_ch, pan_ref_ch).set{index_ch}
  pangenie(reads_ch.combine(ref_ch).combine(index_ch))
}
