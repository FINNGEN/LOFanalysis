version development

workflow regenie_lof {

  input {
    Array[String] phenos
    Float max_maf
    Float info_filter
    String docker
    String prefix
    Int mlogp_filter
    File sets
    File mask
    File lof_variants
    File lof_bgen
    Array[String] lof_list
  }

  scatter(pheno in phenos) {
    call regenie {
      input :
        pheno = pheno,
        prefix = prefix,
        lof_variants = lof_variants,
        sets = sets,
        mask = mask,
        bins = max_maf,
        lof_bgen = lof_bgen,
        mlogp_filter = mlogp_filter
    }
  }

  call merge_results {
    input:
      docker = docker,
      regenie_results = regenie.results,
      comp_files = regenie.sig_res,
      prefix = prefix,
      logs = regenie.log,
      sets = sets,
      max_maf = max_maf,
      info = info_filter,
      lof_list = lof_list,
      mlogp_filter = mlogp_filter
  }

}


task merge_results {

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

  String res_file = prefix + "_lof.txt" # file with all hits
  String sig_file = prefix + "_sig.txt" # file with sig hits (and variant comparison)
  String log_file = prefix + "_lof.log" # file with merged logs
  String sql_file = prefix + "_lof.sql.txt" # file for sql import 
  String readme   = prefix + "_lof_readme" # output readme 
  String var_file = prefix + "_lof_variants.txt" # list of variants used
  String checksum = prefix + "_check.log" # check that all ran or not
  

  command <<<
  # filter results
  paste <(echo PHENO) <(zcat  ~{regenie_results[0]} | head -n2 | tail -n1 | tr ' ' '\t')   > ~{res_file} # write header
  while read f
  do pheno=$(basename $f .regenie.gz |sed 's/~{prefix}_lof_//g' )  && zcat  $f | sed -E 1,2d |  awk -v pheno="$pheno" '{print pheno" "$0}' |  tr ' ' '\t'   >> tmp.txt
  done 	< ~{write_lines(regenie_results)}
  
  sort -grk 13 tmp.txt | grep -vw TEST_FAIL  >> ~{res_file}
  cat tmp.txt | grep -w TEST_FAIL > ~{log_file}
  # merge logs
  while read f
  do paste <( grep -q  "Elapsed" $f && echo 1 || echo 0 ) <(grep "phenoColList"  $f |  awk '{print $2}') <(echo $f| sed 's/\/cromwell_root/gs:\//g')  >> ~{checksum} && cat $f >> ~{log_file}
  done < ~{write_lines(logs)}

  head -n1 ~{comp_files[0]} > ~{sig_file}
  cut -f1,4 ~{sets} > ~{var_file}
  
  while read f
  do

      cat $f | sed -E 1d | awk '{print $5-$6"\t"$0}'   >> sig_tmp.txt
  done <  ~{write_lines(comp_files)}
  sort -rgk 1 sig_tmp.txt | cut -f2- | sed -e 's/ /\t/g' >> ~{sig_file}

  touch ~{sql_file}
  #GZIP FILES
  bgzip ~{log_file} -c > ~{log_file}.gz
  bgzip ~{res_file} -c > ~{res_file}.gz
  
  python3 <<EOF
  import sys,os,gzip,re
  import numpy as np

  # read in gene dict
  with open('~{sets}') as i:
    gene_dict = {}
    for line in i:
        gene,*_,variants = line.strip().split()
        gene_dict[gene] = variants
        
  with open('~{sql_file}','wt') as o, open('~{res_file}') as i:
    next(i) # go past header
    rel = re.findall(r'\d+','~{prefix}')[0] #get release number from prefix
    for line in i:
        pheno,_,_,gene,_,_,_,n,_,beta,se,_,mlogp,_ = line.strip().split()
        try:
            mlogp=float(mlogp)
            if float(mlogp) > float('~{mlogp_filter}'):
                gene = gene.split(".Mask")[0]
                data = map(str,[rel,pheno,gene,gene_dict[gene],np.power(10,-float(mlogp)),beta,se,n])
                o.write('"' + '","'.join(data) + '"\n')  
  
        except:
            print(pheno,gene,mlogp)
            
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
    preemptible: 1
    maxRetries: 2
  }
   
  output {
    File all_hits = res_file +'.gz'
    File sig_hits = sig_file
    File sql_hits = sql_file
    File variants = var_file
    File read     = readme 
    File log      = log_file + '.gz'
    File check    = checksum
  }  
}

task regenie{

  input {
    String docker
    File lof_bgen
    File cov_file
    String pheno
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

  Int mem = ceil(size(lof_bgen,"GB"))*4*cpus
  String regenie_results = prefix + "_lof_" + pheno + ".regenie.gz"

  File sumstats = sub(sumstats_root,"PHENO",pheno)
  File tabix = sumstats + ".tbi"
  String pheno_comp = prefix + "_" + pheno + "_lof_sig_res.txt"

  
  command <<<
  CPUS=$(grep -c ^processor /proc/cpuinfo)

  # REGENIE
  time regenie --step 2 --out ./~{prefix}_lof --threads $CPUS  ~{bargs}  --bgen ~{lof_bgen} --sample ~{lof_sample} --phenoFile ~{cov_file} --covarFile ~{cov_file} --phenoColList ~{pheno} --covarColList ~{covariates} --aaf-bins ~{bins}  --build-mask ~{mask_type} --mask-def ~{mask} --set-list ~{sets} --anno-file ~{lof_variants}
  echo ~{cpus} ~{mem} $CPUS

  ####### SIG HITS FILTER
  
  # KEEP ONLY RELEVANT HITS (MLOGP > FILTER)
  echo "HITS WITH MLOGP > THRESHOLD"
  TMP_HITS="tmp_hits.txt"
  paste <(zcat ~{regenie_results} | awk '$12>~{mlogp_filter}'  | cut -d " " -f 3 | cut -d . -f 1 | awk '{print $0"_GENESTRING"}') <(zcat ~{regenie_results} | awk '$12>~{mlogp_filter}' | cut -d " " -f 1,9,12) | sed -e 's/ /\t/g' > $TMP_HITS

  # GET GENE TO VAR MAPPING FOR ALL GENES WITH HTIS
  cut -f 1,2 ~{lof_variants} | awk '{print $2"_GENESTRING\t"$1}' | sort -k1 | grep -wf <(cut -f1 $TMP_HITS ) > sig_genes.txt
  wc -l sig_genes.txt $TMP_HITS
  
  # FOR SIG GENES GET DATA BY TABIXING FROM SUMSTATS
  echo "RELEVANT GENES WITH SIG MLOGP"
  TMP_VAR="tmp_sig_genes_variants.txt"
  rm -f $TMP_VAR && touch $TMP_VAR
  while read -a  line # for each gene that is a significant HIT tabix from sumstats data
  do
      GENE=${line[0]}
      CHROM=${line[2]}
      POS=${line[3]}
      REF=${line[4]}
      ALT=${line[5]}
      paste <(echo -e $GENE"_GENESTRING\tchr"$CHROM"_"$POS"_"$REF"_"$ALT) <(tabix ~{sumstats} $CHROM:$POS-$POS | grep -w $REF | grep -w $ALT | cut -f 8,9,11) >> $TMP_VAR
  done < <(cat  sig_genes.txt |  sed 's/chr//g' | tr '_' '\t' )

  # THIS IS THE PART TO BE JOINED WITH THE REGENIE RESULT (I COLLAPSE ALL VARIANT DATA FOR EACH GENE INTO A SINGLE LINE)
  TMP_GENE_VAR_HIT="tmp_sig_hits_gene.txt"
  rm -f $TMP_GENE_VAR_HIT && touch $TMP_GENE_VAR_HIT
  while read GENE
  do
      TOP_MLOGP=$( grep -w $GENE $TMP_VAR|  cut -f3 | sort -gr  | head -n1 ) # grep gene and keep highest MLOGP
      GENE_DATA=$(cat $TMP_VAR | grep -w $GENE | cut -f2- | tr '\t' ',' | tr '\n' ','| sed 's/.$//' ) # combine all lines into one
      echo -e $GENE'\t'$TOP_MLOGP'\t'$GENE_DATA >> $TMP_GENE_VAR_HIT
  done < <(cut -f1 $TMP_VAR | sort | uniq)

  # JOIN REGENIE RESULT WITH GENE VARIANT DATA
  PHENO_GENE_HITS="results.txt"
  echo -e "PHENO\tGENE\tCHROM\tGENE_BETA\tGENE_MLOGP\tTOP_VAR_MLOGP\tVAR1,MLOGP1,BETA1,AF1,VAR2..." > $PHENO_GENE_HITS 
  PHENO=~{pheno}
  join -t $'\t' <(sort -k1 $TMP_HITS) <(sort -k1 $TMP_GENE_VAR_HIT) | sed 's/_GENESTRING//g' | awk '{print "PNAME\t"$0}' | sed -e "s/PNAME/${PHENO}/g"  >> $PHENO_GENE_HITS

  mv $PHENO_GENE_HITS ~{pheno_comp}

  >>>

    runtime {
    docker: "~{docker}"
    cpu: "~{cpus}"
    disks:  "local-disk ~{disk_size} HDD"
    memory: "~{mem} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 1
    maxRetries: 2
  }
  
  output {
    File results = regenie_results
    File sig_res = pheno_comp
    File log = prefix + "_lof.log"
  }
}