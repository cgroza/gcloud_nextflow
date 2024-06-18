process check_sickle {
  cpus 1
  memory "2G"
  publishDir "${params.out}", mode: 'copy'

  input:
  tuple val(sample), path(vcf)

  output:
  val(sample), path("${sample}".sickle)

  script:
  """
  bcftools view -H -r chr11:5227001-5227003 ${vcf} > ${sample}.sickle
  """
}

workflow {
  ch_vcfs = Channel.fromPath(params.vcfs).splitCsv(header:true).map{row ->
    [row.sample, file(row.path, checkIfExists:true)]}

  check_sickle(ch_vcfs)
}
