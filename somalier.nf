process extract_bam {
  cpus params.cpus
  memory params.memory
  publishDir params.out

  input:
  tuple path(bam), path(ref), path(sites)

  output:
  path("extracted/*")

  script:
  """
  mkdir extracted
  somalier extract -d extracted/ --sites ${sites} -f ${ref} ${bam}
  """
}

workflow {
  ch_ref = Channel.fromPath(params.ref)
  ch_sites = Channel.fromPath(params.sites)
  ch_bams = Channel.fromPath(params.bams).splitCsv(header:true).map{row -> [file(row.path, checkIfExists:false), file(row.index, checkIfExists:false)]}
  extract_bam(ch_bams.combine(ch_ref).combine(ch_sites))
}
