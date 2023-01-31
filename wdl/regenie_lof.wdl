version development

workflow regenie_lof {

  input {
    String docker
    String prefix
    Float max_maf
    Array[String] chrom_list
  }

  # get lof variants
  call extract_variants { input: docker = docker, max_maf = max_maf}
  # subset vcf to lof variants in each chrom
  scatter (chrom in chrom_list){
    call convert_vcf {input: docker = docker,chrom=chrom,lof_variants = extract_variants.lof_variants }
  }
  
  # merge lof chroms + build vcf/bgen
  call merge    { input: docker = docker, vcfs = convert_vcf.chrom_lof_vcf}

  call create_chunks   { input : docker = docker}
  call build_set_mask  { input : docker = docker, lof_variants=extract_variants.lof_variants}

  scatter ( chunk in create_chunks.all_chunks) {
    call regenie {
      input :
      chunk = chunk,
      prefix=prefix,
      lof_variants = extract_variants.lof_variants,
      sets=build_set_mask.sets,
      mask=build_set_mask.mask,
      bins = max_maf,
      lof_bgen = merge.lof_bgen
    }
  }
   
}


task regenie{

  input {
    String docker
    File lof_bgen
    File cov_file
    File chunk
    File sets
    File mask
    File lof_variants
    String covariates
    String bargs
    String prefix
    Int cpus
    String bins
    String mask_type
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

  time regenie --step 2 --out ./~{prefix}_lof --threads $CPUS  ~{bargs}  --bgen ~{lof_bgen} --sample ~{lof_sample} --pred pred.txt    --phenoFile ~{cov_file} --covarFile ~{cov_file} --phenoColList ~{sep="," phenos} --covarColList ~{covariates} --aaf-bins ~{bins}  --build-mask ~{mask_type} --mask-def ~{mask} --set-list ~{sets} --anno-file ~{lof_variants}

  echo ~{cpus} ~{mem} $CPUS
  wc -l pred.txt 

  >>>

    runtime {
    docker: "~{docker}"
    cpu: "~{cpus}"
    disks:  "local-disk ~{disk_size} HDD"
    memory: "~{mem} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible : 1
  }
  
  output {
    Array[File] results = glob("./*regenie.gz")
    File log = prefix + "_lof.log"
  }
}



task merge {

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
    String docker
    File annot_file
    Array[String] lof_list
    Float max_maf
    Float info_filter
  }

  Int disk_size = ceil(size(annot_file,'GB'))*2 + 10
  Float upper_maf = 1- max_maf
  
  command <<<

  # GET COL NAMES AND NUMBERS
  GIND=$(zcat ~{annot_file} | head -n1 |   tr '\t' '\n' | nl -nln |  grep -w "gene_most_severe" | cut -f1)
  MIND=$(zcat ~{annot_file} | head -n1 |   tr '\t' '\n' | nl -nln |  grep -w "most_severe" | cut -f1)
  IIND=$(zcat ~{annot_file} | head -n1 |   tr '\t' '\n' | nl -nln |  grep -w "INFO" | cut -f1)
  AIND=$(zcat ~{annot_file} | head -n1 |   tr '\t' '\n' | nl -nln |  grep -w "AF" | cut -f1)
  
  echo $GIND $MIND $IIND $AIND

  zcat ~{annot_file} | awk -v OFS='\t' -v c1=$AIND -v c2=$IIND -v c3=$GIND -v c4=$MIND '{print $1,$c1,$c2,$c3,$c4}' | awk '$2 < ~{max_maf} || $2 > ~{1-max_maf}' | awk '$3 > ~{info_filter}' |  grep -wf ~{write_lines(lof_list)} |  tr ':' '_' | awk '{print "chr"$0}' | cut -f 1,4,5 > lof_variants.txt
  
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


task create_chunks{
  input {
    String docker
    File null_map
    Int chunk_phenos
  }
  command <<<
  cat ~{null_map}  >  tmp.txt
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

  output { Array[File] all_chunks = glob("./chunk*") }
}

task build_set_mask {
  
  input {
    String docker
    File lof_variants
  }
  
  command <<<
  
  while read GENE
  do 	GENE_DATA=$( cat ~{lof_variants} | grep -E "(^|[[:space:]])$GENE(\$|[[:space:]])"  | cut -f1 | tr '\n' ',' ) && paste <(echo $GENE) <(echo $GENE_DATA| head -n1 | tr '_' '\t' | sed 's/chr//g' | cut -f -2) <(echo $GENE_DATA| sed 's/.$//' )  >> ./sets.tsv
  done < <(cut -f2 ~{lof_variants} | sort  | uniq )
  
  paste <( echo "Mask1") <(cut -f 3 ~{lof_variants}  | sort | uniq | tr '\n' ',' | sed 's/.$//') > ./mask.txt
  
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
    File mask = "./mask.txt"
    File sets = "./sets.tsv"
  }
}
