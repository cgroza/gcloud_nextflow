params.reads          = "reads.csv"
params.reference      = "reference.fa"
params.vcf            = "variants.vcf"
params.out            = "out"
params.index          = false

process preprocess {
  cpus 20
  memory "120G"

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
  memory "30G"
  publishDir "${params.out}/genotypes"

  input:
  tuple val(sample_name), path(sample_bam), path(cram_ref), path("index/")

  output:
  path("${sample_name}_genotyping.vcf")

  script:
  """
  PanGenie -t 10 -j 10 -s ${sample_name} -i <(samtools fastq --reference ${cram_ref} ${sample_bam}) -f index/index/processed -o ${sample_name}
  """
}

workflow {
  Channel.fromPath(params.cram_reference).set{ref_ch}
  Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path, checkIfExists:false)]}.set{reads_ch}

  if(params.index) {
    Channel.fromPath(params.index).set{index_ch}
  }
  else{
    Channel.fromPath(params.vcf).set{vcf_ch}
    Channel.fromPath(params.pangenome_reference).set{pan_ref_ch}
    preprocess(vcf_ch, pan_ref_ch).set{index_ch}
  }
  pangenie(reads_ch.combine(ref_ch).combine(index_ch))
}
