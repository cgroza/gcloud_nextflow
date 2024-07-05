params.ref = null
params.vcf = null

process make_graph {
  cpus params.cpus
  memory params.memory
  input:
  path(vcf)
  path(fasta)

  output:
  path("index")

  script:
  """
  tabix ${vcf}
  mkdir index
  vg autoindex --tmp-dir \$PWD  -p index/index -w giraffe -v ${vcf} -r ${fasta}
  vg snarls index/index.giraffe.gbz > index/index.pb
  """
}

workflow {
  ch_fasta = Channel.fromPath(params.ref)
  ch_vcf = Channel.fromPath(params.vcf)
  make_graph(ch_vcf, ch_fasta)
}
