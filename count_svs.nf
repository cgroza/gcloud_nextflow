process count_svs {
  cpus 2
  memory "12G"

  input:
  tuple val(sample), path(vcf)

  output:
  path("${sample}.csv")

  script:
  """
  echo -n ${sample}, > ${sample}.csv
  grep -Ev '#|0/0' ${vcf} | wc -l >> ${sample}.csv
  """
}

process merge_counts {
  cpus 1
  memory "6G"
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
