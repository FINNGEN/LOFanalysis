version development

workflow regenie_lof {

  input {
    File null_map
    Array[String] lof_list
    String prefix
    File pheno_file
    File cov_file

    String covariates
    Float max_maf
    Boolean is_binary
  }

  String bio_docker = "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8"
  String regenie_docker = "eu.gcr.io/finngen-refinery-dev/regenie:3.3_r12_cond"

  # read in phenotypes
  Array[String] phenos = transpose(read_tsv(null_map))[0]
  call validate_inputs{input: phenolist=phenos, covariates = covariates, cov = cov_file,pheno = pheno_file, is_binary=is_binary,docker=bio_docker}

  call extract_variants { input: docker = bio_docker, max_maf = max_maf,lof_list=lof_list}
  # subset vcf to lof variants in each chrom
  scatter (chrom in extract_variants.chrom_list){
    call convert_vcf {input: docker = bio_docker,chrom=chrom,lof_variants = extract_variants.lof_variants }
  }

  call merge  { input: docker = bio_docker, vcfs = convert_vcf.chrom_lof_vcf}

  Map[String,File] pheno_null_map =read_map(null_map)
  scatter ( pheno in validate_inputs.validated_phenotypes) {
    call regenie {
      input :
      pheno = pheno,
      docker = regenie_docker,
      cov_file = validate_inputs.validated_cov_pheno_file,
      covariates = validate_inputs.validated_covariates,
      null = pheno_null_map[pheno],
      prefix=prefix,
      lof_variants = extract_variants.lof_variants,
      sets=extract_variants.sets,
      mask=extract_variants.mask,
      bins = max_maf,
      lof_bgen = merge.lof_bgen,
      is_binary = is_binary
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
    String covariates
    File null
    File sets
    File mask
    File lof_variants
    String prefix
    String bins
    Boolean is_binary
  }


  String mask_type = "max"
  String bargs = "--bsize 400 --gz --firth --approx --pThresh 0.01 --firth-se --ref-first"
  
  Int disk_size = ceil(size(lof_bgen,"GB"))*2 + 10
  File lof_index =  lof_bgen + ".bgi"
  File lof_sample = lof_bgen + ".sample"
  Map [String,File] pheno_map = {pheno:null}

  Int cpus = 2
  Int mem = ceil(size(lof_bgen,"GB"))*4*cpus
  String regenie_results = prefix + "_lof_" + pheno + ".regenie.gz"
  String bin_qt = if is_binary == true then "--bt" else "--qt"

  
  command <<<
  CPUS=$(grep -c ^processor /proc/cpuinfo)

  # REGENIE
  cat ~{write_map(pheno_map)} > pred.txt
  time regenie --step 2 ~{bin_qt} --out ./~{prefix}_lof --threads $CPUS  ~{bargs}  --bgen ~{lof_bgen} --sample ~{lof_sample} --pred pred.txt    --phenoFile ~{cov_file} --covarFile ~{cov_file} --phenoColList ~{pheno} --covarColList ~{covariates} --aaf-bins ~{bins}  --build-mask ~{mask_type} --mask-def ~{mask} --set-list ~{sets} --anno-file ~{lof_variants}
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
  
  cut -f2 sets.tsv  | sed 's/chr//g' | sort | uniq | sort -V > chrom_list.txt
  
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
    Array[String] chrom_list = read_lines("chrom_list.txt")

    }
}


task validate_inputs {
  input {
    Array[String] phenolist
    String covariates
    File cov
    File pheno
    Boolean is_binary
    Int minimum_samples_in_group = 5
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
    python3 validate.py ~{cov} ~{pheno} ~{is_binary} "~{sep=","  phenolist}" "~{covariates}" ~{minimum_samples_in_group}

    # join files
    PHENOCOLS=$(cat pheno_columns_to_take)
    join --header -1 1 -2 1 -t $'\t' \
        <(cat <(zcat -f ~{cov}|head -n1) <(zcat -f ~{cov}|tail -n+2|sort -k1)) \
        <(cat <(zcat -f ~{pheno}|head -n1) <(zcat -f ~{pheno}|tail -n+2|sort -k1)|cut -f "$PHENOCOLS") \
gzip > validated_cov_pheno.tsv.gz 

    
  >>>
  runtime {
    preemptible: 2
    disks: "local-disk 20 HDD"
    docker: "${docker}"
    cpu: 1
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "4 GB"
  }
