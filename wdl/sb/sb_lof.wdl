version development

workflow regenie_lof {

  input {
    File null_map
    Array[String] lof_list
    String prefix
    File cov_file
    Float max_maf
  }

  String bio_docker = "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8"
  String regenie_docker = "eu.gcr.io/finngen-refinery-dev/regenie:3.3_r12_cond"

  call validate_inputs {input :null_map = null_map,cov_file = cov_file,prefix = prefix,lof_list = lof_list, docker = bio_docker}

  call extract_variants { input: docker = bio_docker, max_maf = max_maf,lof_list=lof_list}

  # subset vcf to lof variants in each chrom
Array[String] chrom_list =  ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
  scatter (chrom in chrom_list){
    call convert_vcf {input: docker = bio_docker,chrom=chrom,lof_variants = extract_variants.lof_variants }
  }

  call merge  { input: docker = bio_docker, vcfs = convert_vcf.chrom_lof_vcf}

  scatter ( p_data in read_tsv(null_map)) {
    String pheno = p_data[0]
    call regenie {
      input :
      docker = regenie_docker,
      cov_file = cov_file,
      pheno = pheno,
      null = p_data[1],
      prefix=prefix,
      lof_variants = extract_variants.lof_variants,
      sets=extract_variants.sets,
      mask=extract_variants.mask,
      bins = max_maf,
      lof_bgen = merge.lof_bgen,
    }
  }

  call merge_results {
    input:
    docker = bio_docker,
    regenie_results = regenie.results,
    prefix=prefix,
    logs = regenie.log,
    sets=extract_variants.sets,
  }

  output {
    File all_hits = merge_results.all_hits
    File var_file = merge_results.variants
    File log_file = merge_results.log
    File check    = merge_results.check
  }

}


task merge_results {

  input {
    String docker
    Array[File] regenie_results
    Array[File] logs
    File sets
    # README STUFF
    String prefix
  }

  String res_file = prefix + "_lof.txt" # file with all hits
  String log_file = prefix + "_lof.log" # file with merged logs
  String var_file = prefix + "_lof_variants.txt" # list of variants use
  String checksum = prefix + "_check.log" # check that all ran or not
  
  command <<<
  # filter results
  paste <(echo PHENO) <(zcat  ~{regenie_results[0]} | head -n2 | tail -n1 | tr ' ' '\t')   > ~{res_file} # write header
  while read f
  do pheno=$(basename $f .regenie.gz |sed 's/~{prefix}_lof_//g' )  && zcat  $f | sed -E 1,2d |  awk -v pheno="$pheno" '{print pheno" "$0}' |  tr ' ' '\t'   >> tmp.txt
  done 	< ~{write_lines(regenie_results)}
  
  sort -grk 13 tmp.txt | grep -vw TEST_FAIL >> ~{res_file}
  cat tmp.txt | grep -w TEST_FAIL > ~{log_file}
  # merge logs
  while read f
  do paste <( grep -q  "Elapsed" $f && echo 1 || echo 0 ) <(grep "phenoColList"  $f |  awk '{print $2}') <(echo $f| sed 's/\/cromwell_root/gs:\//g')  >> ~{checksum} && cat $f >> ~{log_file}
  done < ~{write_lines(logs)}

  # variants used
  cut -f1,4 ~{sets} > ~{var_file}
  >>>
 runtime {
    docker: "~{docker}"
    cpu: 1
    disks:  "local-disk 10 HDD"
    memory: "2 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible : 1
  }
   
  output {
    File all_hits = res_file
    File log      = log_file
    File check    = checksum
    File variants = var_file
    
  }
}
    
 
 
task regenie{

  input {
    String docker
    File lof_bgen
    File cov_file
    String pheno
    File null
    File sets
    File mask
    File lof_variants
    String prefix
    String bins
  }


  String covariates = "SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC{1:10},IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm_relift,BATCH_DS18_MIGRAINE_1_norm_relift,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm_relift,BATCH_DS21_SUPER_2_norm_relift,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm_nosymmetric,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm_relift,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm"

  String mask_type = "max"
  String bargs = "--bt --bsize 400 --gz --firth --approx --pThresh 0.01 --firth-se --ref-first"
  
  Int disk_size = ceil(size(lof_bgen,"GB"))*2 + 10
  File lof_index =  lof_bgen + ".bgi"
  File lof_sample = lof_bgen + ".sample"
  Map [String,File] pheno_map = {pheno:null}

  Int cpus = 2
  Int mem = ceil(size(lof_bgen,"GB"))*4*cpus
  String regenie_results = prefix + "_lof_" + pheno + ".regenie.gz"


  
  command <<<
  CPUS=$(grep -c ^processor /proc/cpuinfo)

  # REGENIE
  cat ~{write_map(pheno_map)} > pred.txt
  time regenie --step 2 --out ./~{prefix}_lof --threads $CPUS  ~{bargs}  --bgen ~{lof_bgen} --sample ~{lof_sample} --pred pred.txt    --phenoFile ~{cov_file} --covarFile ~{cov_file} --phenoColList ~{pheno} --covarColList ~{covariates} --aaf-bins ~{bins}  --build-mask ~{mask_type} --mask-def ~{mask} --set-list ~{sets} --anno-file ~{lof_variants}
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
    File results = regenie_results
    File log = prefix + "_lof.log"
  }
}




task merge {

  input {
    Array[File] vcfs
    String docker
  }

  String bargs =  "-filetype vcf -bgen-bits 8 -bgen-compression zlib -vcf-genotype-field GP -bgen-permitted-input-rounding-error 0.005 -ofiletype bgen_v1.2 "
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
    File lof_variants
    String chrom
    String docker
  }

  String vcf_root = "gs://finngen-production-library-red/finngen_R12/genotype_2.0/data/finngen_R12_chrCHROM.vcf.gz"
  File vcf = sub(vcf_root,"CHROM",chrom)
  File index = vcf + ".tbi"
  Int disk_size = ceil(size(vcf,"GB"))*2
  
  String out_file = chrom + "_lof.vcf.gz"
  String log_file = chrom + "_lof.log"
  
  command <<<
  # create list of positions and list of variants fixing chrom 23/X issues. Variants are labelled X but chrom file is under 23
  cut -f1 ~{lof_variants} | sed 's/chrX/chr23/g' | grep chr~{chrom}  |  sed 's/chr23/chrX/g' > chrom_lof_variants.txt
  head chrom_lof_variants.txt
  cat chrom_lof_variants.txt | awk -F "_" '{print $1"\t"$2"\t"$2}' > chrom_lof_positions.txt
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
    Array[String] lof_list
    Float max_maf
    Float info_filter
    Int gene_variants_min_count
  }

  File annot_file = "gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/annotations/R12_annotated_variants_v0.small.gz"
  Int disk_size = ceil(size(annot_file,'GB'))*2 + 10
  Float upper_maf = 1- max_maf
  
  command <<<
  # GET COL NAMES AND NUMBERS
  GIND=$(zcat ~{annot_file} | head -n1 |   tr '\t' '\n' | nl -nln |  grep -w "gene_most_severe" | cut -f1)
  MIND=$(zcat ~{annot_file} | head -n1 |   tr '\t' '\n' | nl -nln |  grep -w "most_severe" | cut -f1)
  IIND=$(zcat ~{annot_file} | head -n1 |   tr '\t' '\n' | nl -nln |  grep -w "INFO" | cut -f1)
  AIND=$(zcat ~{annot_file} | head -n1 |   tr '\t' '\n' | nl -nln |  grep -w "AF" | cut -f1)
  
  echo $GIND $MIND $IIND $AIND

  zcat ~{annot_file} | awk -v OFS='\t' -v c1=$AIND -v c2=$IIND -v c3=$GIND -v c4=$MIND '{print $1,$c1,$c2,$c3,$c4}' | awk '$2 < ~{max_maf} || $2 > ~{1-max_maf}' | awk '$3 > ~{info_filter}' |  grep -wf ~{write_lines(lof_list)} |  tr ':' '_' | awk '{print "chr"$0}' | cut -f 1,4,5 > tmp.txt

  # keep only genes with >1 variants
  cat tmp.txt| awk '{print $1"\t"$2"_GENESTRING\t"$3}' | grep -wf <(cut -f2 tmp.txt | sort | uniq -c | awk '{$1=$1;print}' | awk '$1>=~{gene_variants_min_count}' | cut -d " " -f2 | sort -k1 | awk '{print $1"_GENESTRING"}') | sed 's/_GENESTRING//g'  > lof_variants.txt
  
  while read GENE
  do 	GENE_DATA=$( cat lof_variants.txt | grep -E "(^|[[:space:]])$GENE(\$|[[:space:]])"  | cut -f1 | tr '\n' ',' ) && paste <(echo $GENE) <(echo $GENE_DATA| head -n1 | tr '_' '\t' | sed 's/chr//g' | cut -f -2) <(echo $GENE_DATA| sed 's/.$//' )  >> ./sets.tsv
  done < <(cut -f2 lof_variants.txt | sort  | uniq )
  
  paste <( echo "Mask1") <(cut -f 3 lof_variants.txt  | sort | uniq | tr '\n' ',' | sed 's/.$//') > ./mask.txt
  

  
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
    File mask = "./mask.txt"
    File sets = "./sets.tsv"

    }
}


task validate_inputs {

  input {
    File null_map
    Array[String] lof_list
    String prefix
    File cov_file
    String docker
  }

  Int disk = ceil(size(cov_file,'GB'))*2
  Int mem = disk + 10


  command <<<

  # GREP ALL STRINGS INTO A SINGLE FILE AND THEN CHECK FOR ERROR THERE
  grep -vP '^[a-zA-Z0-9_.-]*$' <( cut -f1 ~{null_map}) >> tmp.txt
  grep -vP '^[a-zA-Z0-9_.-]*$' ~{write_lines(lof_list)} >> tmp.txt
  grep -vP '^[a-zA-Z0-9_.-]*$' <( echo ~{prefix}) >> tmp.txt
  grep -vP '^[a-zA-Z0-9_.-]*$' <(zcat -f ~{cov_file} | head -n1 | tr "\t" "\n") >> tmp.txt

  set -o pipefail
  if grep -q . tmp.txt; then
      echo "Irregular pattern found " >&2
      cat tmp.txt >&2
      exit 1
  else
      :
  fi

  # write all columns to file
  cut -f1 ~{null_map} > phenos.txt
  # python code that checks if columns are valid in terms of case count
  python3 <<CODE
  import sys
  import pandas as pd
  import numpy as np

  pheno_file='~{cov_file}'
  with open("phenos.txt") as i:cols=[elem.strip() for elem in i.readlines()]
  print(f"{len(cols)} columns left to analyze")

  # read first 10k lines to filter out most problematic phenos and float
  df = pd.read_csv(pheno_file,engine='python',usecols=cols,sep=None,index_col=None,nrows = 1000).select_dtypes(include='number').astype(bool)

  # cols that do not pass the test
  check_cols = list(df.columns[np.where(df.sum(axis=0).lt(5))])
  print(f"{len(check_cols)} columns left to analyze")
  if not check_cols:sys.exit(0) #exit without error if all good
  
  df = pd.read_csv(pheno_file,sep="\t",index_col=None,usecols=check_cols,dtype=float).astype(bool)
  cols = list(df.columns[np.where(df.sum(axis=0).lt(5))])
  if cols:print(f"{cols} do not have at least 5 cases")

  #invert logically DF and check for controls
  or_cols = list(df.columns[np.where((~df).sum(axis=0).lt(5))])
  if or_cols:print(f"{or_cols} do not have at least 5 controls")

  print(len(cols),len(or_cols))
  if cols+or_cols:sys.exit("Error")

  CODE
  >>>

    runtime {
    docker: "~{docker}"
    cpu: 1
    disks:  "local-disk ~{disk} HDD"
    memory: "~{mem} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible : 1
  }
}
