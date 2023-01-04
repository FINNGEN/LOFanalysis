version 1.0

workflow filter_lof{

  input {
    String docker
    Array[String] chrom_list
    String bargs
  }

  # get lof variants
  call extract_variants { input: docker = docker}

  # subset vcf to lof variants in each chrom
  scatter (chrom in chrom_list){
    call convert_vcf {input: docker = docker,chrom=chrom,lof_variants = extract_variants.lof_variants }
  }

  # merge lof chroms + build vcf/bgen
  call merge_vcf    { input: docker = docker, vcfs = convert_vcf.chrom_lof_vcf,bargs=bargs}
  
}


task merge_vcf {

  input {
    Array[File] vcfs
    String docker
    String bargs
  }
  
  Int disk_size = ceil(size(vcfs,'GB'))*2 + 20
  String out_file = "lof"
  
  command <<<
  # VCF CONCATENATION
  cat ~{write_lines(vcfs)} | sort -g  > vcf_list.txt 
  bcftools concat -f vcf_list.txt -Oz -o ~{out_file}.vcf.gz
  tabix ~{out_file}.vcf.gz
  bcftools index -n ~{out_file}.vcf.gz
  # BGEN CONVERSION
  qctool -g ~{out_file}.vcf.gz -og ~{out_file}.bgen  -os ~{out_file}.bgen.sample  ~{bargs}
  bgenix -g ~{out_file}.bgen -clobber -index
  >>>

  runtime {
    docker: "~{docker}"
    cpu: "4"
    disks:   "local-disk ~{disk_size} HDD"
    memory: "16 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    
  }
    
  output{
    File lof_bgen   = out_file + ".bgen"
    File lof_sample = out_file + ".bgen.sample"
    File lof_index  = out_file + ".bgen.bgi"
    File lof_vcf    = out_file + ".vcf.gz"
    File lof_vcf_tbi= out_file + ".vcf.gz.tbi"
    
  }

}
  
task convert_vcf {

  input {
    String vcf_root
    File lof_variants
    String chrom
    String docker
  }

  File vcf = sub(vcf_root,"CHROM",chrom)
  File index = vcf + ".tbi"
  Int disk_size = ceil(size(vcf,"GB"))*2
  
  String out_file = chrom + "_lof.vcf.gz"
  String log_file = chrom + "_lof.log"
  
  command <<<
  # create list of positions and list of variants fixing chrom 23/X issues. Variants are labelled X but chrom file is under 23
  cut -f1 ~{lof_variants} | sed 's/chrX/chr23/g' | grep chr~{chrom}  |  sed 's/chr23/chrX/g' > chrom_lof_variants.txt
  head chrom_lof_variants.txt
  cat chrom_lof_variants.txt |  awk -F "_" '{print $1"\t"$2"\t"$2}' > chrom_lof_positions.txt
  head chrom_lof_positions.txt
  
  # create new vcf
  echo "building vcf ..."
  bcftools view ~{vcf} -R chrom_lof_positions.txt --i ID=@chrom_lof_variants.txt  -Oz -o ~{out_file}

  # build index
  echo "building index ..."
  tabix ~{out_file}
  # sanity check that all variants are there
  paste <(wc -l < chrom_lof_variants.txt) <(bcftools index -s ~{out_file} | cut -f 3) > ~{log_file}

  >>>
  
  runtime {
    docker: "~{docker}"
    cpu: "4"
    disks:   "local-disk ~{disk_size} HDD"
    memory: "16 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
  }
  
  output {
    File chrom_lof_vcf = out_file
    File chrom_lof_vcf_index  = out_file + ".tbi"
    File chrom_log = log_file
    } 

}


task extract_variants {

  input {
    File annotated_variants
    Array[String] lof_list
    String docker

    File freq_file
    Float max_maf
    
    File info_score_file
    String info_filter
  }

  Int disk_size = ceil(size(info_score_file,"GB") + size(freq_file,"GB")) + 10
  Float upper_maf = 1- max_maf
  
  command <<<

  #get info score
  zcat ~{info_score_file} | sed -E 1d |  cut -f1,2 | awk '$2 >  ~{info_filter}' | cut -f1 | sort > info_filter.txt
  wc -l info_filter.txt
  # MAF filter
  cat ~{freq_file} | sed -E 1d | cut -f 2,5 | awk '$2 < ~{max_maf} || $2 > ~{upper_maf}' | cut -f1 | sort > maf_filter.txt
  wc -l maf_filter.txt
  # shared variants
  comm -12 maf_filter.txt info_filter.txt > pass_variants.txt
  wc -l pass_variants.txt
  # lof variants
  zcat ~{annotated_variants} | cut -f 3,5,6 | grep -wf ~{write_lines(lof_list)} | sort -k1  > lof_tmp.txt
  join -t $'\t' pass_variants.txt lof_tmp.txt > lof_variants.txt

  >>>
  runtime {
    docker: "~{docker}"
    cpu: "1"
    disks:   "local-disk ~{disk_size} HDD"
    memory: "2 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
  }
  
  output {
    File lof_variants = "lof_variants.txt"
    }

}

