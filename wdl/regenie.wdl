workflow LOF_regenie{

    Array[String] chrom_list
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

    call split_phenos {}

    scatter (pheno_list in split_phenos.pheno_lists){
      call step2 {
        input:
        null_list = pheno_list,
        n_phenos = split_phenos.n_phenos,
        pgen = merge_pgen.pgen,
        psam = merge_pgen.psam,
        pvar = merge_pgen.pvar
      }

    }
}


task step2{

    File null_list
    Array[File] nulls = read_lines(null_list)
    String out_root = "/cromwell_root/finngen_lof" + basename(null_list)

    File covarFile
    String covariates
    File pgen
    File psam
    File pvar

    File mask
    File sets
    File annotation

    String docker
    Int n_phenos
    Int disk_size =  2*ceil(size(pgen,'GB'))+ n_phenos*(ceil(size(nulls[0],"GB"))) +10


    command <<<
    find /cromwell_root/ -name '*loco.gz'  > local_null_list.txt
    while read f; do basename  -s .loco.gz $f >> phenos.txt ; done < local_null_list.txt
    paste phenos.txt local_null_list.txt > pred_list.txt
    cat pred_list.txt

    pheno_cols=`cut -f1 pred_list.txt |paste -s -d,`

    regenie \
    --step 2 \
    --bt  \
    --pgen ${sub(pgen,"lof.pgen","lof")} \
    --ref-first \
    --covarFile ${covarFile} \
    --covarColList ${covariates} \
    --pred pred_list.txt \
    --phenoFile  ${covarFile} \
    --phenoColList $pheno_cols \
    --set-list ${sets} \
    --anno-file ${annotation} \
    --mask-def ${mask} \
    --threads ${n_phenos} \
    --firth \
    --approx \
    --aaf-bins 0.01,0.1,0.5 \
    --bsize 200 \
    --out ${out_root}

    >>>

    output {
    Array[File] results = glob("/cromwell_root/*regenie")
    File run_log = out_root +".log"

    }
    runtime {
        docker: "${docker}"
        cpu: n_phenos
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
