version development

workflow regenie_lof {

  input {
    File null_map
    Float max_maf
    Float info_filter
    Array[String] lof_list
    Array[String] chrom_list
    String docker
    String prefix
    Int mlogp_filter
    
  }

  # get lof variants
  call extract_variants { input: docker = docker, max_maf = max_maf,info_filter=info_filter,lof_list=lof_list}
  # subset vcf to lof variants in each chrom
  scatter (chrom in chrom_list){
    call convert_vcf {input: docker = docker,chrom=chrom,lof_variants = extract_variants.lof_variants }
  }
  
  # merge lof chroms + build vcf/bgen
  call merge    { input: docker = docker, vcfs = convert_vcf.chrom_lof_vcf}



  Array[Array[String]] pheno_data = read_tsv(null_map)
  
  scatter ( p_data in pheno_data) {
    String pheno = p_data[0]
    call regenie {
      input :
      pheno = pheno,
      null = p_data[1],
      prefix=prefix,
      lof_variants = extract_variants.lof_variants,
      sets=extract_variants.sets,
      mask=extract_variants.mask,
      bins = max_maf,
      lof_bgen = merge.lof_bgen,
      mlogp_filter = mlogp_filter
    }
  }

  call merge_results {input: docker = docker,regenie_results = regenie.results,comp_files = regenie.sig_res,prefix=prefix,logs = regenie.log,sets=extract_variants.sets,max_maf=max_maf,info=info_filter,lof_list=lof_list,mlogp_filter=mlogp_filter}

}


task merge_results{

  input {
    String docker
    Array[File] regenie_results
    Array[File] comp_files
    Array[File] logs
    File sets
    # README STUFF
    String prefix
    File lof_template
    String max_maf
    Int mlogp_filter
    Array[String] lof_list
    String info

  }

  String res_file = prefix + "_lof.txt"
  String sig_file = prefix + "_sig.txt"
  String log_file = prefix + "_lof.log"
  String sql_file = prefix + "_lof.txt.gz"
  String readme   = prefix + "_lof_readme"
  String var_file = prefix + "_lof_variants.txt"
  String checksum = prefix + "_check.log"
  

  command <<<
  # filter results
  paste <(echo PHENO) <(zcat  ~{regenie_results[0]} | head -n2 | tail -n1 | tr ' ' '\t')   > ~{res_file} # write header
  while read f
  do pheno=$(basename $f .regenie.gz |sed 's/~{prefix}_lof_//g' )  && zcat  $f | sed -E 1,2d |  awk -v pheno="$pheno" '{print pheno" "$0}' |  tr ' ' '\t'   >> tmp.txt
  done 	< ~{write_lines(regenie_results)}

  sort -grk 13 tmp.txt >> ~{res_file}
  # merge logs
  while read f
  do paste <( grep -q  "Elapsed" $f && echo 1 || echo 0 ) <(grep "phenoColList"  $f |  awk '{print $2}') <(echo $f| sed 's/\/cromwell_root/gs:\//g')  >> ~{checksum} && cat $f >> ~{log_file}
  done < ~{write_lines(logs)}

  head -n1 ~{comp_files[0]} > ~{sig_file}

  cut -f1,4 ~{sets} > ~{var_file}
  
  while read f
  do

      cat $f | sed -E 1d | awk '{print $5-$6"\t"$0}' >> sig_tmp.txt
  done <  ~{write_lines(comp_files)}
  sort -rgk 1 sig_tmp.txt | cut -f2- >> ~{sig_file}
  
  python3 <<EOF
  import sys,os,gzip
  import numpy as np
  sets='~{sets}'
  hits='~{res_file}'
  
  # read in gene dict
  with open(sets) as i:
    gene_dict = {}
    for line in i:
        gene,*_,variants = line.strip().split()
        gene_dict[gene] = variants
  
  
  out_file = '~{sql_file}'
  with gzip.open(out_file,'wt') as o, open(hits) as i:
    next(i)
    o.write('\t'.join(['PHENO','GENE','variants','p.value','BETA','SE','N']) +'\n')
    for line in i:
        pheno,_,_,gene,_,_,_,n,_,beta,se,_,mlogp,_ = line.strip().split()
        pval = str(np.power(10,-float(mlogp)))
        gene = gene.split(".Mask")[0]
        variants = gene_dict[gene]
        out_line = '\t'.join([pheno,gene,variants,pval,beta,se,n]) + '\n'
        o.write(out_line)
    

  lof_list = f"[{'~{sep="," lof_list}'}]"
  tags = [("[PREFIX]",'~{prefix}'),("[N_GENES]",str(len(gene_dict))),("[MAF]",'~{max_maf}'),("[INFO]",'~{info}'),('[LOF_LIST]',lof_list),('[MLOGP]','~{mlogp_filter}')]
  with open('~{lof_template}') as i,open('~{readme}','wt') as o:
    for line in i:
        for tag in tags:
            line = line.replace(tag[0],tag[1])
        o.write(line)

  EOF  

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
    File sig_hits = sig_file
    File sql_hits = sql_file
    File variants = var_file
    File read     = readme 
    File log      = log_file
    File check    = checksum
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
    String covariates
    String bargs
    String prefix
    Int cpus
    String bins
    String mask_type
    String sumstats_root
    Int mlogp_filter
  }

  Int disk_size = ceil(size(lof_bgen,"GB"))*2 + 10
  File lof_index =  lof_bgen + ".bgi"
  File lof_sample = lof_bgen + ".sample"
  Map [String,File] pheno_map = {pheno:null}

  Int mem = ceil(size(lof_bgen,"GB"))*4*cpus
  String regenie_results = prefix + "_lof_" + pheno + ".regenie.gz"

  File sumstats = sub(sumstats_root,"PHENO",pheno)
  File tabix = sumstats + ".tbi"
  String pheno_comp = prefix + "_" + pheno + "_lof_sig_res.txt"

  
  command <<<
  CPUS=$(grep -c ^processor /proc/cpuinfo)
  cat ~{write_map(pheno_map)} > pred.txt

  time regenie --step 2 --out ./~{prefix}_lof --threads $CPUS  ~{bargs}  --bgen ~{lof_bgen} --sample ~{lof_sample} --pred pred.txt    --phenoFile ~{cov_file} --covarFile ~{cov_file} --phenoColList ~{pheno} --covarColList ~{covariates} --aaf-bins ~{bins}  --build-mask ~{mask_type} --mask-def ~{mask} --set-list ~{sets} --anno-file ~{lof_variants}

  echo ~{cpus} ~{mem} $CPUS
  wc -l pred.txt

  # KEEP ONLY GENES >1 VARIANT COUNT
  echo "SIG GENES"
  cut -f2 ~{lof_variants} | sort | uniq -c | awk '{$1=$1;print}' | awk '$1>1' | cut -d " " -f2 | sort -k1 | join -2 2 - <(cut -f 1,2 ~{lof_variants}| sort -k 2 ) | awk '{print $1"_GENESTRING\t"$2}'> tmp_genes.txt
  wc -l tmp_genes.txt
  
  # KEEP ONLY RELEVANT HITS (LOW PVAL)
  echo "HITS WITH MLGOP > 6"
  TMP_HITS="tmp_hits.txt"
  paste <(zcat ~{regenie_results} | awk '$12 >6' | cut -d " " -f 3 | cut -d . -f 1 | awk '{print $0"_GENESTRING"}') <(zcat ~{regenie_results} | awk '$12>~{mlogp_filter}' | cut -d " " -f 1,9,12) > $TMP_HITS
  wc -l $TMP_HITS
  
  # GET INTERSECTION OF SIG HITS AND RELEVANT GENES, TABIXING DATA FROM SUMSTATS
  echo "RELEVANT GENES WITH SIG MLOGP"
  TMP=tmp_variants.txt
  rm -f $TMP && touch $TMP
  while read -a line
  do
      CHROM=${line[0]}
      POS=${line[1]}
      ALT=${line[3]}
      tabix ~{sumstats} $CHROM:$POS-$POS | grep -w $ALT | cut -f 8,9,11 >> $TMP
  done < <(cat tmp_genes.txt | grep -wf <(cut -f1 $TMP_HITS ) | cut -f2 | sed 's/chr//g' | tr '_' '\t' )

  # add gene and variant id to file above
  TMP_VAR="tmp_sig_genes_variants.txt"
  paste <(cat tmp_genes.txt | grep -wf <(cut -f1 $TMP_HITS) ) $TMP > $TMP_VAR
  wc -l $TMP_VAR
  
  # HERE I BUILD THE ACTUAL LINE OF TEXT THAT NEEDS TO BE PASTED,MUNGING THE ABOVE FILE
  TMP_GENE_VAR_HIT="tmp_sig_hits_gene.txt"
  rm -f $TMP_GENE_VAR_HIT && touch $TMP_GENE_VAR_HIT
  while read GENE
  do
      echo $GENE
      TOP_MLOGP=$(sort -rk3 $TMP_VAR | grep -w $GENE |cut -f3 | head -n1 )
      GENE_DATA=$(cat $TMP_VAR | grep -w $GENE | cut -f2- | tr '\t' ',' | tr '\n' ','| sed 's/.$//' )
      echo -e $GENE'\t'$TOP_MLOGP'\t'$GENE_DATA >> $TMP_GENE_VAR_HIT
  done < <(cut -f1 $TMP_VAR | sort | uniq)

  # JOIN REGENIE RESULT WITH VARIANT DATA
  PHENO_GENE_HITS="results.txt"
  echo -e "PHENO\tGENE\tCHROM\tGENE_BETA\tGENE_MLOGP\tTOP_VAR_MLOGP\tVAR1,MLGOP1,BETA1,AF1,VAR2..." > $PHENO_GENE_HITS 

  PHENO=~{pheno}
  join <(sort -k1 $TMP_HITS) <(sort -k1 $TMP_GENE_VAR_HIT) | sed 's/_GENESTRING//g' | awk '{print "PNAME\t"$0}' | sed -e "s/PNAME/${PHENO}/g"  >> $PHENO_GENE_HITS

  mv $PHENO_GENE_HITS ~{pheno_comp}

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
    File sig_res = pheno_comp
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

  zcat ~{annot_file} | awk -v OFS='\t' -v c1=$AIND -v c2=$IIND -v c3=$GIND -v c4=$MIND '{print $1,$c1,$c2,$c3,$c4}' | awk '$2 < ~{max_maf} || $2 > ~{1-max_maf}' | awk '$3 > ~{info_filter}' |  grep -wf ~{write_lines(lof_list)} |  tr ':' '_' | awk '{print "chr"$0}' | cut -f 1,4,5 > tmp.txt

  # keep only genes with >1 variants
  cat tmp.txt| awk '{print $1"\t"$2"_GENESTRING\t"$3}' | grep -wf <(cut -f2 tmp.txt | sort | uniq -c | awk '{$1=$1;print}' | awk '$1>1' | cut -d " " -f2 | sort -k1 | awk '{print $1"_GENESTRING"}') | sed 's/_GENESTRING//g'  > lof_variants.txt
  
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



