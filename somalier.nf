process ancestry {
  cpus params.cpus
  memory params.memory
  publishDir params.out

  input:
  path(samples)
  path(queries)
  path(labels)

  output:
  path("somalier-ancestry*")

  script:
  """
  somalier ancestry --labels ${labels} ${samples} ++ ${queries}
  """
}
process extract_bam {
  cpus params.cpus
  memory params.memory
  publishDir params.out

  input:
  tuple path(bam), path(bam_index), path(ref), path(sites)

  output:
  path("extracted/*")

  script:
  """
  mkdir extracted
  somalier extract -d extracted/ --sites ${sites} -f ${ref} ${bam}
  """
}

workflow {
  if(params.extract) {
    ch_ref = Channel.fromPath(params.ref)
    ch_sites = Channel.fromPath(params.sites)
    ch_bams = Channel.fromPath(params.bams).splitCsv(header:true).map{row -> [file(row.path, checkIfExists:false), file(row.index, checkIfExists:false)]}
    extract_bam(ch_bams.combine(ch_ref).combine(ch_sites))
  }
  if(params.ancestry) {
    ch_labels = Channel.fromPath(params.labels)
    ch_samples = Channel.fromPath(params.samples).collect()
    ch_queries = Channel.fromPath(params.queries).collect()
  }
}
