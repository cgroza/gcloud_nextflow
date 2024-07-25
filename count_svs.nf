process count_svs {
  cpus 1
  memory "2G"

  input:
  tuple val(sample), path(vcf)

  output:
  path("${sample}.csv")

  script:
  """
  echo -n ${sample}, > ${sample}.csv
  grep -v '0/0' ${vcf} >> ${sample}.csv
  """
}

process merge_counts {
  cpus 1
  memory "2G"
  publishDir "${params.out}", mode: 'copy'

  input:
  path(counts)

  output:
  path("counts.csv")

  script:
  """
  cat ${counts} > counts.csv
  """
}

workflow {
  ch_vcfs = Channel.fromPath(params.vcfs).splitCsv(header:true).map{row ->
    [row.sample, file(row.path, checkIfExists:true)]}

  merge_counts(count_svs(ch_vcfs).collect())
}
