workflow LOF_regenie{

    Array[String] chrom_list
    Boolean test
    scatter (chrom in chrom_list){
        call convert_vcf_pgen{
            input:
                chrom = chrom}

        }

    call sum {input: values = convert_vcf_pgen.size}
    call merge_pgen {
        input:
            sizes =  sum.out,
            pgens = convert_vcf_pgen.pgen,
            pvars = convert_vcf_pgen.pvar,
            psams = convert_vcf_pgen.psam
            }

    call split_phenos {input: test = test}

    scatter (pheno_list in split_phenos.pheno_lists){
      call step2 {
        input:
        null_list = pheno_list,
        chrom_list = chrom_list,
        test = test,
        pgen = merge_pgen.pgen,
        psam = merge_pgen.psam,
        pvar = merge_pgen.pvar
      }

    }
}


task step2{

    File null_list
    Array[File] nulls = read_lines(null_list)
    String out_root = "/cromwell_root/finngen_lof/"
    Boolean test

    Array[String] chrom_list
    Array[String] chrom = if test then ["22"] else chrom_list
    String regenie_args
    Int max_retries
    String name

    File covarFile
    String covariates
    File pgen
    File psam
    File pvar

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
     --chrom ${sep=',' chrom}  --regenie-args ${regenie_args}  \
     --max-retries ${max_retries}  --name ${name}

    >>>

    output {
    Array[File] results = glob("${out_root}/*regenie")
    Array[File] logs = glob("${out_root}/*log")
    File failed = out_root + name + "_failed.txt"
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpus}"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        memory: "8 GB"
    }

}
task split_phenos{

    File null_list
    Int split_size

    String docker
    Boolean test
    String test_filter = if test then "| head -n " + split_size*2 else ""

    command <<<
    cat ${null_list} ${test_filter} | split -l ${split_size} -d - pheno_list
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
    }
}

task merge_pgen {

    Array[File] pgens
    Array[File] pvars
    Array[File] psams

    Int sizes
    Int disk_size = sizes*3 + 10

    String docker

    command <<<

    while read f; do echo $f | cut -f 1 -d '.'  >> merge_list.txt ; done < <(cat ${write_lines(pgens)} |  sort -V)
    cat merge_list.txt
    plink2 --pmerge-list merge_list.txt --out lof
    >>>
    runtime {
        docker: "${docker}"
        cpu: 16
        disks: "local-disk " + "${disk_size}" + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        memory: "8 GB"
    }
    output {
        File pgen = "./lof.pgen"
        File pvar = "./lof.pvar"
        File psam = "./lof.psam"
       }
}

task convert_vcf_pgen {

    String chrom_root
    String chrom

    File chrom_vcf = sub(chrom_root,'CHROM',chrom)
    File chrom_tbi = chrom_vcf + ".tbi"
    File lof_variants
    String plink_params

    String disk_size = 2*ceil(size(chrom_vcf,'GB'))+10
    String docker

    command <<<
    plink2 --vcf ${chrom_vcf}  --extract ${lof_variants}  ${plink_params}  --make-pgen  --out ${chrom}_lof
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


task sum {

  Array[Float] values
  String docker

  command <<<
    python -c "print(int(${sep="+" values}))"
  >>>
  output {
    Int out = read_int(stdout())
  }
  runtime {
    docker: "${docker}"
    memory: "1 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    disks: "local-disk 10 HDD"
    cpu: 1
  }
}
