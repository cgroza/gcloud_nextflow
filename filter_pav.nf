process filter_svs {
  cpus 1
  memory "2G"
  publishDir "${params.out}", mode: 'copy'

  input:
  tuple val(sample), path(vcf)

  output:
  tuple val(sample), path("${sample}.sickle")

  script:
  """
  bcftools view -i 'INFO/SVLEN>=50 || INFO/SVLEN<=-50' -Oz -o ${sample}.vcf.gz ${vcf} 
  """
}

workflow {
  ch_vcfs = Channel.fromPath(params.vcfs).splitCsv(header:true).map{row ->
    [row.sample, file(row.path, checkIfExists:true)]}

  filter_svs(ch_vcfs)
}
