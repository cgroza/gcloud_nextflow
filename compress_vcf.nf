process compress_vcf {
  cpus params.cpus
  memory params.memory

  input:
  path(vcf)

  output:
  tuple path("${vcf}.gz"), path("${vcf}.gz.tbi")

  script:
  """
  bgzip ${vcf}
  tabix ${vcf}.gz
  """
}

workflow {
  ch_vcfs = compress_vcf(Channel.fromPath(params.vcf))
}
