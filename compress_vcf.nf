process compress_vcf {
  cpus params.cpus
  memory params.memory
  publishDir params.out, mode: "copy"

  input:
  path(vcf)

  output:
  tuple path("${vcf}.gz"), path("${vcf}.gz.tbi")

  script:
  """
  bcftools sort ${vcf}| bgzip --threads ${params.cpus} > ${vcf}.gz
  tabix ${vcf}.gz
  """
}

workflow {
  ch_vcfs = compress_vcf(Channel.fromPath(params.vcf))
}
