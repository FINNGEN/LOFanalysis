version 1.0

workflow gp_lof{

  input {
    File lof_vcf
    String docker
    File lof_map
  }

  call build_gp_vcf    { input : docker = docker, lof_map=lof_map,lof_vcf = lof_vcf }
  call bgen_conversion { input : docker = docker, lof_vcf = build_gp_vcf.gene_vcf}
}


task bgen_conversion{

  input {
    File lof_vcf
    String  docker
  }

  File lof_vcf_index = lof_vcf + ".tbi"
  Int disk_size = ceil(size(lof_vcf,"GB"))*10 + 20
  String out_root = basename(lof_vcf,".vcf.gz")
  
  command <<<
  qctool -g ~{lof_vcf} -og ~{out_root}.bgen -os ~{out_root}.sample -vcf-genotype-field GP -filetype vcf
  bgenix -g ~{out_root}.bgen -clobber -index
  >>>
  
  runtime {
    docker: "~{docker}"
    cpu: 8
    disks:   "local-disk ~{disk_size} HDD"
    memory: "2 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
  }
  output {
    File lof_bgen   = out_root + ".bgen"
    File lof_index  = out_root + ".bgen.bgi"
    File lof_sample = out_root + ".sample"
  }
}


task build_gp_vcf {

  input {
    File lof_vcf
    File lof_map
    String docker
    Int cpus   
  }
  
  File lof_index = lof_vcf + ".tbi"
  Int disk_size = ceil(size(lof_vcf,"GB"))*10 + 20
  
  command <<<
  df -h .
  python3.7 /scripts/lof_gp.py --vcf ~{lof_vcf} --lof_map ~{lof_map} -o .
  df -h .
  >>>

  runtime {
    docker: "~{docker}"
    cpu:    "~{cpus}"
    disks:  "local-disk ~{disk_size} HDD"
    memory: "2 GB"
    zones:  "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
  }

  output {
    File gene_vcf   = "./lof_gene_gp.vcf.gz"
    File gene_index = "./lof_gene_gp.vcf.gz.tbi"
    File gene_map   = "./lof_gene_gp_dict.tsv"
  }

}
