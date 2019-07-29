workflow LOF_saige{

     File chrom_info
     String docker
     Array[Array[String]] chrom_info_list = read_tsv(chrom_info)
     String lof
     File samples


     call filter_variants{
     	  input:
		docker = docker,
		lof = lof
          }
     scatter (data in chrom_info_list){

     	     call gene_matrix{
	     	  input:
		  lof = lof,
		  chrom_data = data,
		  docker = docker,
		  samples = samples,
		  lof_variants = filter_variants.variant_gene_list
	     }
     }

     call merge_matrix {

     	  input:
	  lof = lof,
	  matrix_files = gene_matrix.matrix_file,
	  docker =docker
     }
     
     File variance_list
     call return_phenotypes {
     	  input:
     	  variance_list = variance_list,
	  docker = docker
	  }
	  
     scatter (pheno in return_phenotypes.phenotypes){
     	     call pheno_saige {
	     	  input:
		  matrix = merge_matrix.gene_matrix,
		  pheno = pheno,
		  samples = samples
	     	     }
     	     }

     }


task filter_variants{

     String lof_plink_root 
     File lof_bed = lof_plink_root  + ".bed"
     File lof_bim = lof_plink_root  + ".bim"
     File lof_fam = lof_plink_root  + ".fam"

     File lof_variant_genes

     String lof
     Float max_maf
     Float geno
     String docker

     String dollar = "$"

     command {
     plink_root=$( echo ${lof_bed} | sed 's/.bed//g' )
     plink2 --bfile ${dollar}plink_root --max-maf ${max_maf} --geno ${geno} --write-snplist --out ${lof}_variants
     join <(sort ${lof_variant_genes}) <(sort ${lof}_variants.snplist) -t ${dollar}'\t' > ${lof}_gene_list.txt

     }

     output {
     File variant_gene_list = "${lof}_gene_list.txt"
     }
     runtime {
        docker: "${docker}"
        cpu: 4
        memory: "1 GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    	}        

}

task merge_results{

     Array[File] gene_results
     File header
     String docker
     
     command <<<
     cat ${header} > gene_results.txt &&\
     while read f ; do  awk '{gsub(".*/","",FILENAME);gsub(".SAIGE.txt","",FILENAME);print FILENAME" "$0}' $f |  sed -E 1d  >> tmp.txt; done < ${write_lines(gene_results)} && sort -gk 9 tmp.txt| awk '$9>0' >> gene_results.txt && bgzip gene_results.txt
     
     >>>
     output {
     File results = "gene_results.txt.gz"
     }
     runtime {
        docker: "${docker}"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    	}
     }


task return_phenotypes{

     File variance_list
     String docker
        
     command <<<
     cat ${variance_list} | xargs basename -s .varianceRatio.txt | awk -F- '{ print $NF}' > phenotypes.txt
     >>>
     output {	
	Array[String] phenotypes = read_lines("phenotypes.txt")
	}
	
     runtime {
        docker: "${docker}"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 1 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    	}
}

task pheno_saige {
	
	File matrix
	
	String pheno
	String variance_basepath 
	File variance_file = sub(variance_basepath,"PHENO",pheno) 
	File null_file= sub(variance_file, ".varianceRatio.txt", ".rda")
	File samples
	
	Int minmac
	String loco

	Int mem
	String docker
	Int cpu

	String outfile =  pheno +  ".SAIGE.txt"
	
	
	command {
	step2_SPAtests.R \
	--dosageFile=${matrix} \
	--dosageFileNrowSkip=1 \
	--dosageFileNcolSkip=1 \
	--dosageFilecolnamesSkip="GENE" \
	--minMAC=${minmac} \
	--sampleFile=${samples} \
	--LOCO=${loco} \
	--numLinesOutput=1000 \
	--IsOutputAFinCaseCtrl=TRUE \
	--GMMATmodelFile=${null_file} \
	--varianceRatioFile=${variance_file} \
	--SAIGEOutputFile=${outfile} 
	}
	
	output {	
	File out = outfile
	}
	
	runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    	}
}




task gene_matrix {

     Int disk_factor
     String chrom_file
     Array[String] chrom_data
     
     String chrom = chrom_data[0]
     File chrom_vcf = sub(chrom_file,'CHROM',chrom)
     File chrom_tbi = chrom_vcf + ".tbi"
     File lof_variants
     File samples

     # get chrom size and convert to disk size
     Int chrom_size = chrom_data[1]
     String disk_size =disk_factor * chrom_size

     Float call_filter
     String call_type = if call_filter >0 then "--hard-call " + call_filter else "--gp0"
     
     Float? info_score
     String info_filter = if defined(info_score) then "--info_score " + info_score else " "

     String lof          
     String docker
     Int cpu

     command{
     python3 /Scripts/lof.py \
     --lof_variants ${lof_variants} \
     --vcf ${chrom_vcf} \
     -o /cromwell_root/ \
     -c ${chrom} \
     ${info_filter} \
     --lof ${lof} \
     ${call_type} \
     -s ${samples} \
     --cpus ${cpu}
     }

     runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        zones: "europe-west1-b"
        preemptible: 1
	bootDiskSizeGb: 20
	memory: "8 GB"
    }
    output {	
       File matrix_file = "/cromwell_root/${chrom}/${lof}_${chrom}_gene_matrix.tsv"
       File gene_dict = "/cromwell_root/${chrom}/${lof}_${chrom}_gene_variants_dict.txt"
       }

}

task merge_matrix {

     Array[File] matrix_files
     String docker
     String name
     String lof
     
     command <<<
     ls /cromwell_root/ &&
     awk 'FNR==1 && NR!=1 { while (/^GENE/) getline; } 1 {print}' ${sep=" " matrix_files} > /cromwell_root/${name}_${lof}_gene_matrix.tsv
     >>>
     
     runtime {
        docker: "${docker}"
        cpu: "1"
	disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
        preemptible: 1
	bootDiskSizeGb: 20
	memory: "4 GB"
    }

    output{
    File gene_matrix = "/cromwell_root/${name}_${lof}_gene_matrix.tsv"
    }
}