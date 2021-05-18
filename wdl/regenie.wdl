

workflow LOF_regenie{

    Array[String] chrom_list
    Boolean test
    Float max_maf
    scatter (chrom in chrom_list){
        call convert_vcf_pgen{
            input:
              max_maf = max_maf,
              chrom = chrom}

        }

    call merge_pgen {
        input:
            pgens = convert_vcf_pgen.pgen,
            pvars = convert_vcf_pgen.pvar,
            psams = convert_vcf_pgen.psam
            }

    call split_phenos {input: test = test}

    call annotate {
      input :
        freq = merge_pgen.freq
    }

    Array[File] pheno_lists = split_phenos.pheno_lists
    String name
    scatter (idx in range(length(pheno_lists))) {
      call step2 {
        input:
        annotation = annotate.annot,
        sets = annotate.sets,
        null_list = pheno_lists[idx],
        max_maf = max_maf,
        chrom_list = chrom_list,
        test = test,
        name = name + "_" + idx,
        pgen = merge_pgen.pgen,
        psam = merge_pgen.psam,
        pvar = merge_pgen.pvar,
        freq = merge_pgen.freq
      }

    }
    call merge_results{
      input:
        name = name,
        max_maf = max_maf,
        regenie_results = step2.results,
        header = step2.header[0],
        regenie_logs = step2.log
    }
}

task merge_results{


  Array[File] regenie_results
  Array[File] regenie_logs
  Float max_maf
  File header
  String docker
  String name


  command <<<
  cat ${header} | bgzip -c > ${name}_results.txt.gz
  while read f; do cat $f | sed 's/.most_severe.${max_maf}//g'   >> tmp.txt ; done < ${write_lines(regenie_results)}

  cat tmp.txt |sort -rgk 12 | bgzip -c >> ${name}_results.txt.gz

  while read f; do cat $f   >> ${name}_logs.txt ; done < ${write_lines(regenie_logs)}
  >>>
  output {
    File merged_log = "./${name}_logs.txt"
    File merged_results = "./${name}_results.txt.gz"
  }
  runtime {
      docker: "${docker}"
      cpu: 1
      disks: "local-disk 10 HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      memory: "2 GB"
      preemptible: 2
  }
}




task step2{

    File null_list
    Array[File] nulls = read_lines(null_list)
    String out_root = "/cromwell_root/finngen_lof/"
    Boolean test

    Array[String] chrom_list
    Float max_maf
    String regenie_args
    Int max_retries
    String name

    File covarFile
    String covariates
    File pgen
    File psam
    File pvar
    File freq

    File mask
    File sets
    File annotation

    String docker
    Int disk_size =  2*ceil(size(pgen,'GB'))+ length(nulls)*(ceil(size(nulls[0],"GB"))) +10
    Int cpus

    command <<<
    find /cromwell_root/ -name '*loco.gz'  > local_null_list.txt
    while read f; do basename  -s .loco.gz $f >> phenos.txt ; done < local_null_list.txt
    paste phenos.txt local_null_list.txt > pred_list.txt
    cat pred_list.txt

    python3 /Scripts/regenie.py -o ${out_root} --covariates ${covariates} \
     --annot ${annotation} --pred pred_list.txt  --pgen ${pgen} \
     --mask ${mask}  --sets ${sets}  --pheno-file ${covarFile}  \
     --chrom ${sep=' ' chrom_list}  --regenie-args ${regenie_args}  \
     --max-retries ${max_retries} \
     --aaf-bins ${max_maf} --aaf-file ${freq} \
     --name ${name} ${if test then "--test "  else ""}

    >>>

    output {
    File log = out_root + name + ".log"
    File results = out_root + name + ".regenie"
    File header =  out_root + name + ".header"
    File failed  = out_root + name + "_failed.txt"
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpus}"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        memory: "8 GB"
    }

}


task annotate{

  File annotation
  File freq
  String docker
  Int disk_size = ceil(size(annotation,"GB")) + 10

  command <<<
  python3 /Scripts/annotate_regenie.py -o /cromwell_root/ --annot ${annotation} --variant-file ${freq}

  >>>

  output {
    File annot = "./annot.tsv"
    File sets = "./sets.tsv"

  }

  runtime {
      docker: "${docker}"
      cpu: 1
      disks: "local-disk ${disk_size} HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      memory: "2 GB"
      preemptible: 2
  }

}

task split_phenos{

    File null_list
    Int split_size

    String docker
    Boolean test

    command <<<
    cat ${null_list} ${if test then "| head -n " + split_size*2 else ""} | split -l ${split_size} -d - pheno_list
    >>>
    output {
      Array[File] pheno_lists = glob("/cromwell_root/pheno_list*")
      Int n_phenos = split_size
    }
    runtime {
        docker: "${docker}"
        cpu: 1
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        memory: "8 GB"
        preemptible: 2
    }
}

task merge_pgen {

    Array[File] pgens
    Array[File] pvars
    Array[File] psams

    Int size = ceil(size(pgens[10],'GB')) * length(pgens) * 3 + 10
    String docker

    command <<<

    while read f; do echo $f | cut -f 1 -d '.'  >> merge_list.txt ; done < <(cat ${write_lines(pgens)} |  sort -V)
    cat merge_list.txt
    plink2 --pmerge-list merge_list.txt --out lof
    plink2 --pfile lof --freq --out lof
    cut -f 2,5 lof.afreq | sed -E 1d > lof.freq
    >>>
    runtime {
        docker: "${docker}"
        cpu: 16
        disks: "local-disk " + "${size}" + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        memory: "8 GB"
        preemptible: 2
    }
    output {
        File pgen = "./lof.pgen"
        File pvar = "./lof.pvar"
        File psam = "./lof.psam"
        File freq = "./lof.freq"
       }
}

task convert_vcf_pgen {

    String chrom_root
    String chrom

    File chrom_vcf = sub(chrom_root,'CHROM',chrom)
    File chrom_tbi = chrom_vcf + ".tbi"
    File lof_variants
    String plink_params
    Float max_maf

    String disk_size = 2*ceil(size(chrom_vcf,'GB'))+10
    String docker

    command <<<
    plink2 --vcf ${chrom_vcf}  --extract ${lof_variants}  ${plink_params}  --make-pgen  --out ${chrom}_lof --max-maf ${max_maf}
    >>>

    runtime {
        docker: "${docker}"
        cpu: 4
        disks: "local-disk " + "${disk_size}" + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 1
        memory: "8 GB"
    }
    output {
        File pgen = "${chrom}_lof.pgen"
        File pvar = "${chrom}_lof.pvar"
        File psam = "${chrom}_lof.psam"
        Int size = ceil(size(chrom_vcf,'GB'))
       }
}
