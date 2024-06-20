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
  ls | grep '.vcf.gz' > vcfs.txt
  bcftools merge --threads ${params.merge_threads} -m none -l vcfs.txt -Oz -o merged_bcftools.vcf.gz
  tabix merged_bcftools.vcf.gz
  truvari collapse -i merged_bcftools.vcf.gz -o merged.vcf
  """
}

workflow {
  ch_vcfs = Channel.fromPath(params.vcfs).splitCsv(header:true).map{row ->
    [file(row.path, checkIfExists:true)]}.collect()

  trouvari_merge(ch_vcfs)
}
