process index_vcf {
  cpus 1
  memory "2G"

  input:
  path(vcf)

  output:
  tuple path(vcf), path("${vcf}.tbi")

  script:
  """
  tabix ${vcf}
  """
}

process trouvari_merge {
  cpus params.merge_threads
  memory params.merge_memory
  publishDir "${params.out}", mode: 'copy'

  input:
  path(vcfs)

  output:
  path("merged.vcf"), emit: sv_variants_ch

  script:
  """
  ls *.vcf.gz > vcfs.txt
  bcftools merge --threads ${params.merge_threads} -m none -l vcfs.txt -Oz -o merged_bcftools.vcf.gz
  tabix merged_bcftools.vcf.gz
  truvari collapse -i merged_bcftools.vcf.gz -o merged.vcf
  """
}

workflow {
  ch_vcfs = index_vcf(Channel.fromPath(params.vcfs).splitCsv(header:true).map{row ->
    [file(row.path, checkIfExists:true)]})

  trouvari_merge(ch_vcfs.collect())
}
