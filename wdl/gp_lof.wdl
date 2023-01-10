version development

workflow gp_lof{

  input {
    File lof_vcf
    String docker
    File lof_map
    File null_map
    Boolean test
    Int regenie_cpus
  }

  call build_gp_vcf    { input : docker = docker, lof_map=lof_map,lof_vcf = lof_vcf,test=test}
  call bgen_conversion { input : docker = docker, lof_vcf = build_gp_vcf.gene_vcf}
  call create_chunks   { input : docker = docker, null_map = null_map,test=test,cpus=regenie_cpus}

  scatter ( chunk in create_chunks.all_chunks) {
    call regenie_gp_lof {  input : chunk = chunk, lof_bgen = bgen_conversion.lof_bgen,cpus=regenie_cpus}
  }
}

task regenie_gp_lof {

  input {
    String regenie_docker
    File lof_bgen
    File cov_file
    File chunk
    String covariates
    String bargs
    String prefix
    Int cpus
  }

  Int disk_size = ceil(size(lof_bgen,"GB"))*2 + 10
  
  File lof_index =  lof_bgen + ".bgi"
  File lof_sample = lof_bgen + ".sample"
  
  # this localizes the files
  Map [String,File] pheno_map = read_map(chunk)
  # get list of phenos only
  Array[String] phenos = keys(pheno_map)
  Int mem = ceil(size(lof_bgen,"GB"))*4*cpus
  command <<<
  
  CPUS=$(grep -c ^processor /proc/cpuinfo)

  cat ~{write_map(pheno_map)} > pred.txt
  time regenie --step 2 --out ./~{prefix}_lof --threads $CPUS  ~{bargs}  --bgen ~{lof_bgen} --sample ~{lof_sample} --pred pred.txt    --phenoFile ~{cov_file} --covarFile ~{cov_file} --phenoColList ~{sep="," phenos} --covarColList ~{covariates}

  echo ~{cpus} ~{mem} $CPUS
  wc -l pred.txt 
  >>>
  runtime {
    docker: "~{regenie_docker}"
    cpu: "~{cpus}"
    disks:  "local-disk ~{disk_size} HDD"
    memory: "~{mem} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible : 1
  }
  
  output {
    Array[File] results = glob("./*regenie.gz")
    File logs = prefix + "_lof.log"
  }
}

task build_gp_vcf {

  input {
    File lof_vcf
    File lof_map
    String docker
    Int cpus
    Boolean test
  }

  File lof_index = lof_vcf + ".tbi"
  Int disk_size = ceil(size(lof_vcf,"GB"))*10 + 20
  
  command <<<
  python3 /scripts/lof_gp.py --vcf ~{lof_vcf} --lof_map ~{lof_map} -o . --log info ~{if test then "--test" else ""}
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



task bgen_conversion{

  input {
    File lof_vcf
    String  docker
  }

  File lof_vcf_index = lof_vcf + ".tbi"
  Int disk_size = ceil(size(lof_vcf,"GB"))*10 + 20
  String out_root = basename(lof_vcf,".vcf.gz")
  
  command <<<
  qctool -g ~{lof_vcf} -og ~{out_root}.bgen -os ~{out_root}.bgen.sample -vcf-genotype-field GP -filetype vcf
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
    File lof_sample = out_root + ".bgen.sample"
  }
}

task create_chunks{
  input {
    String docker
    File null_map
    Int chunk_phenos
    Int cpus
    Boolean test
  }


  command <<<
  cat ~{null_map} | shuf  ~{if test then " | head -n " + chunk_phenos*cpus else ""} >  tmp.txt
  split -del ~{chunk_phenos} tmp.txt chunk --additional-suffix=.txt
  >>>
  
  runtime {
    docker: "~{docker}"
    cpu: 2
    disks:  "local-disk 1 HDD"
    memory: "2 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible : 1
  }

  output {
    Array[File] all_chunks = glob("./chunk*")
  }
}



