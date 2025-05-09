params.reads          = "reads.csv"
params.index_dir      = false
params.prefix         = false
params.out            = "out"
params.cpus           = 22
params.memory         = "88G"
params.Q              = 20


process giraffe {
  cpus params.cpus
  memory params.memory
  machineType "c3-standard-22"
  publishDir "${params.out}/packs"

  input:
  tuple val(sample_name), path(sample_bam), path(cram_ref), path(index)

  output:
  tuple path("${sample_name}.gz"), path("${sample_name}.pack")

  script:
  """
  samtools fastq -N -@ ${params.cpus} --reference ${cram_ref} -1 ${sample_name}_1.fq.gz -2 ${sample_name}_2.fq.gz -0 /dev/null -s /dev/null ${sample_bam}
  vg giraffe -t ${params.cpus} -N ${sample_name} --index-basename ${index}/${params.prefix} -f ${sample_name}_1.fq.gz -f ${sample_name}_2.fq.gz | \
    vg pack -d -t ${params.cpus} -Q ${params.Q} -x ${index}/${params.prefix}.gbz -g - -o ${sample_name}.pack | \
    pigz > ${sample_name}.gz
  rm ${sample_name}_1.fq.gz ${sample_name}_2.fq.gz ${cram_ref}.fai
  """
}

process gfa2bin {
  cpus params.cpus
  memory params.memory
  publishDir "${params.out}/plink/"

  input:
  path(packs)

  output:
  path("plink*")

  script:
  """
  ls ${packs} | parallel gunzip
  ls ${packs} | awk -v OFS='\t' '{print(\$1,\$1)}' > packlist
  gfa2bin cov -p packlist -o plink
  """

}

workflow {
  Channel.fromPath(params.cram_reference).set{ref_ch}
  Channel.fromPath(params.reads).splitCsv(header:true).map{row -> [row.sample, file(row.path, checkIfExists:false)]}.set{reads_ch}

  Channel.fromPath("${params.index_dir}", type : "dir", checkIfExists : true).set{index_ch}
  gfa2bin(giraffe(reads_ch.combine(ref_ch).combine(index_ch)).map{s -> p[0]}.collect())
}
