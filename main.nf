process survivor_merge {
  cpus params.merge_threads
  memory params.merge_memory
  publishDir "${params.out}", mode: 'copy'

  input:
  path(vcfs)

  output:
  path("merged.vcf"), emit: sv_variants_ch
  path("vcfs.txt")

  script:
  """
  ls *.vcf > vcfs.txt
  SURVIVOR merge vcfs.txt 0.1 0 1 0 0 100 merged.vcf
  """
}

workflow {
  ch_vcfs = Channel.fromPath(params.vcfs).splitCsv(header:true).map{row ->
    [file(row.path, checkIfExists:true)]}

  survivor_merge(ch_vcfs)
}
