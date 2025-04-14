version development

workflow regenie_lof_exome_custom {

  input {
    String prefix
    File pheno_list

    # OPTIONAL PARAMETERS for building MASK
    File? custom_variants
    File? custom_sets
    File? custom_mask
    # REGENIE PARAMS
    File pheno_file
    File cov_file
    String covariates
    Boolean is_binary
    
  }

  String bio_docker = "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8"
  String regenie_docker = "eu.gcr.io/finngen-refinery-dev/regenie:3.3_r12_cond"

  # all custom inputs into one
  Array[File?] custom_inputs = [custom_variants,custom_sets,custom_mask]
  # if even just one is not defined, we resort back to default
  if (length(select_all(custom_inputs)) != 3) {
    call extract_variants { input: docker = bio_docker}
  
  }
  File lof_variants = select_first([extract_variants.lof_variants,custom_variants])
  File sets = select_first([extract_variants.sets,custom_sets])
  File mask = select_first([extract_variants.mask,custom_mask])
  

  # Array[String] chrom_list =  ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
  Array[String] chrom_list =  ["22"]
  scatter (chrom in chrom_list){
    call convert_vcf {input:docker = bio_docker,chrom=chrom,lof_variants = lof_variants }
  }
  call merge  { input: docker = bio_docker, vcfs = convert_vcf.chrom_lof_vcf}

  
  Array[String] phenos = read_lines(pheno_list)
  call validate_inputs{
    input:
    phenolist=phenos,
    covariates = covariates,
    cov = cov_file,
    pheno = pheno_file,
    is_binary=is_binary,
    docker=bio_docker
  }

  #STEP1

  #STEP2
}




task extract_variants {

  input {
    String docker
    File annot_file
    Array[String] lof_list
    Float max_maf
    Int gene_variants_min_count
    
  }
    
  command <<<
  # Extract variables
  annot_file=~{annot_file}	       
  max_maf=~{max_maf}
  lof_list=~{write_lines(lof_list)}
  gene_variants_min_count=~{gene_variants_min_count}
    
  # Get column indices
  GIND=$(zcat -f "${annot_file}" | head -1 | awk -F'\t' '{for(i=1; i<=NF; i++) if($i == "gene_most_severe") {print i; exit;}}')
  MIND=$(zcat -f "${annot_file}" | head -1 | awk -F'\t' '{for(i=1; i<=NF; i++) if($i == "most_severe") {print i; exit;}}')
  AIND=$(zcat -f "${annot_file}" | head -1 | awk -F'\t' '{for(i=1; i<=NF; i++) if($i == "AF") {print i; exit;}}')
  
  echo "$GIND $MIND  $AIND"
  #SUBSET ONLY TO VARIANTS WITH MAX MAF < THRESHOLD AND WITH LOF VARIANTS
  zcat -f "${annot_file}" | awk -v OFS='\t' -v c1="$AIND"  -v c2="$GIND" -v c3="$MIND" '{print $1,$c1,$c2,$c3}' |   awk -v max_maf="${max_maf}" '$2 > 0 && $2 < max_maf || $2 > 1-max_maf && $2 < 1'|  grep -wf ${lof_list} |  cut -f 1,3,4 |  sort > tmp.txt


  # keep only genes with >1 variants
  awk -F'\t' '{gene=$2; variant=$1; if(!(gene in variants)){variants[gene]=variant; count[gene]=1} else {variants[gene]=variants[gene] "," variant; count[gene]++}} END {for(gene in variants){print gene "\t" count[gene] "\t" variants[gene]}}' tmp.txt | sort | awk -v min_count="${gene_variants_min_count}" '$2>=min_count' > sets.tsv

  # NOW I NEED TO SUBSET THE VARIANTS AND GENERATE THE MASK
  cat sets.tsv  | cut -f3 | tr ',' '\n' | sort > lof_variants.txt
  paste <( echo "Mask1") <(join -t $'\t'  lof_variants.txt  tmp.txt | cut -f 3 | sort | uniq | tr '\n' ',' |  sed 's/,$/\n/') > ./mask.txt
  
  >>>
  runtime {
    docker: "~{docker}"
    cpu: "1"
    disks:   "local-disk ~{ceil(size(annot_file,'GB'))*2 + 10} HDD"
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

task convert_vcf {

  input {
    File lof_variants
    String chrom
    String docker
    String vcf_root
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


task validate_inputs {

  input {
    Array[String] phenolist
    String covariates
    File cov
    File pheno
    Boolean is_binary
    String docker
  }
  output {
    Array[String] validated_phenotypes = read_lines("validated_phenotypes")
    File validated_cov_pheno_file = "validated_cov_pheno.tsv.gz"
    String validated_covariates = read_string("validated_covariates")
  }
  command <<<

  set -eux
  # check phenolist. It should contain only endpoints.
  cat << "__EOF__" > validate.py
  ##
  # This script does for the following:
  # - Check that all covariates are in covariate file
  # - check that all endpoints are in endpoint file
  # - Check that each endpoint has at least minimum_samples_in_group samples, and at least 5 samples, even after joining with covariates.
  #   - Missing values in endpoints OR covariates are filtered out. Missing values are coded as NA.
  # - Take phenotype column indices, for later joining with covariate file. 
  # 
  ##
  import sys
  import re
  import gzip
  from collections import defaultdict as dd
  from contextlib import contextmanager

  @contextmanager
  def uopen(fname,oper_type):
      """
      Universal opener
      Open both gzipped and plaintext files with no fuzz
      """
      gz_magicnumber=b"\x1f\x8b"
      type="normal"
      with open(fname,"rb") as f:
          if f.read()[0:2] == gz_magicnumber:
              type="gz"
      if type=="normal":
          with open(fname,oper_type) as f:
              yield f
      elif type == "gz":
          with gzip.open(fname,oper_type) as f:
              yield f
      else:
          raise Exception("invalid file format")
  ## input values
  cov_file = sys.argv[1]
  pheno_file = sys.argv[2]
  is_binary = sys.argv[3] == "true"
  phenos = sys.argv[4].split(",")
  covariates = sys.argv[5].split(",")
  minimum_samples_in_group = max(int(sys.argv[6]),5)
  ## filter out expandable covariates (e.g. "PC{1:10}"), to be expanded later
  filtered_covars = [a for a in covariates if "{" not in a]
  ## create expanded covariates, and add them to list of variants
  expand_covars = [a for a in covariates if "{" in a]
  for v in expand_covars:
      v_base = v.split("{")[0]
      v_r_start,v_r_end = v.split("{")[1].split("}")[0].split(":")[0:2]
      start = int(v_r_start)
      end = int(v_r_end)
      for i in range(start,end+1):
          filtered_covars.append(f"{v_base}{i}")
  ## check that all phenotypes are alphanumeric+"-_"
  phenomatch = lambda x:bool(re.fullmatch("[a-zA-Z0-9_-]*",x))
  if not all([phenomatch(a) for a in phenos]):
      raise Exception(f"Not all endpoints were alphanumeric! Here are the invalid phenotypes:{[a for a in phenos if not phenomatch(a)]}")
  ## Check that each phenotype has at least minimum_samples_in_group samples, in both cases and controls or in total if quantitative endpoints
  ## Also take in phenotype column numbers now for easier join later
  pheno_colnumbers = []
  with uopen(pheno_file,"rt") as pheno_f, uopen(cov_file,"rt") as cov_f:
      pheno_h = pheno_f.readline().strip().split("\t")
      if not all([a in pheno_h for a in phenos]):raise Exception(f"phenotypes {[a for a in phenos if a not in pheno_h]} not in phenotype file!")
      phdi = {a:i for i,a in enumerate(pheno_h)}
      pheno_colnumbers = [phdi[a]+1 for a in phenos]
      cov_h = cov_f.readline().strip().split("\t")
      # check that each covariate is in cov file
      if not all([a in cov_h for a in filtered_covars]):
          raise Exception(f"covariates {[a for a in filtered_covars if a not in cov_h]} not in covariate file!")
      chdi = {a:i for i,a in enumerate(cov_h)}
      pheno_samples = dd(lambda : dd(set))
      if is_binary:
          for line in pheno_f:
              cols = line.strip().split("\t")
              sample = cols[phdi["IID"]]
              for p in phenos:
                  val = cols[phdi[p]]
                  pheno_samples[p][val].add(sample)
      else:
          for line in pheno_f:
              cols = line.strip().split("\t")
              sample = cols[phdi["IID"]]
              for p in phenos:
                  if cols[phdi[p]] !="NA":
                      pheno_samples[p]["quant"].add(sample)
      #read all samples in covariate file
      #Check that covariates are not NA. Regenie reads "NA" as a missing value, even though it is not being said in the documentation.
      cov_samples = set()
      for line in cov_f:
          add=True
          cols = line.strip().split("\t")
          for c in filtered_covars:
              if cols[chdi[c]] == "NA":
                  add=False
          if add:
              cov_samples.add(cols[chdi["IID"]])
      # for each pheno, check that there are enough samples per endpoint
      bin_groups = {"0":"controls","1":"cases"}
      for p in phenos:
          if is_binary:
              for group in ["0","1"]:
                  n_samples = len(pheno_samples[p][group].intersection(cov_samples))
                  if n_samples < minimum_samples_in_group:
                      raise Exception(f"There are not enough samples with with covariate data for phenotype {p} and group '{bin_groups[group]}'! There are {len(cov_samples)} samples with covariate data, {len(pheno_samples[p][group])} samples for {p} {bin_groups[group]}, but their intersection is only {n_samples} samples, less than minimum allowed {minimum_samples_in_group}")
          else:
              n_pheno_samples = len(pheno_samples[p]["quant"])
              n_samples = len(pheno_samples[p]["quant"].intersection(cov_samples))
              if n_samples < minimum_samples_in_group:
                  raise Exception(f"There are not enough samples with covariate data for quantitative phenotype {p}! There are {len(cov_samples)} samples with covaria tedata, {n_pheno_samples} samples for phenotype {p}, but their intersection is only {n_samples} samples, less than minimum allowed {minimum_samples_in_group}")
  with open("validated_phenotypes","w") as f:
      for p in phenos:
          f.write(f"{p}\n")
  with open("validated_covariates","w") as f:
      f.write(",".join(covariates))
  with open("pheno_columns_to_take","w") as f:
      f.write(",".join(map(str,[1]+pheno_colnumbers)))
  __EOF__
  # validate the shape and types of the phenotype file
  python3 validate.py ~{cov} ~{pheno} ~{is_binary} "~{sep=","  phenolist}" "~{covariates}" 5

  # join files
  PHENOCOLS=$(cat pheno_columns_to_take)
  join --header -1 1 -2 1 -t $'\t' <(cat <(zcat -f ~{cov}|head -n1) <(zcat -f ~{cov}|tail -n+2|sort -k1))  <(cat <(zcat -f ~{pheno}|head -n1) <(zcat -f ~{pheno}|tail -n+2|sort -k1)|cut -f "$PHENOCOLS") | gzip > validated_cov_pheno.tsv.gz 
      
  >>>
  runtime {
    preemptible: 2
    disks: "local-disk 20 HDD"
    docker: "~{docker}"
    cpu: 1
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "4 GB"
  }

}
