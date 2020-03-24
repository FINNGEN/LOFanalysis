 workflow LOF_saige{
    
    Array[String] chrom_list
    String docker
    String lof
    String prefix
    File samples
    File gene_list   
    String variance_list

    # test mode
    Boolean test
    String var_list = if test then sub(variance_list,".txt","_test.txt") else variance_list
     
    # RETURN LIST OF PHENOTYPES TO RUN
    call return_phenotypes {
        input:
	variance_list = var_list,
	docker = docker
    }

    ##############
    #---MATRIX---#
    ##############
    # RETURN GENE MATRIX FOR EACH CHROM
    scatter (chrom in chrom_list){
     	call gene_matrix{
	    input:
	    lof = lof,
	    chrom = chrom,
	    docker = docker,
	    samples = samples,
	    lof_variants = gene_list
	}
    }

    # MERGE MATRICES INTO ONE AND CREATE FINAL GENE TO VARIANT JSON
    call merge_matrix {
     	input:
	lof = lof,
        dicts = gene_matrix.gene_dict,
	vcfs = gene_matrix.vcf,
	docker =docker,
        prefix = prefix
    }

    #############
    #---SAIGE---#
    #############
    # RUN SAIGE FOR EACH PHENO
    scatter (pheno in return_phenotypes.phenotypes){
     	call pheno_saige {
            input:
            gene_vcf = merge_matrix.gene_vcf,
            gene_tbi = merge_matrix.gene_tbi,
            pheno = pheno,
            samples = samples,
        }
    }

    # MERGE SAIGE RESULTS
    call merge_saige{
     	input:
        gene_results = pheno_saige.out,
        docker = docker
    }


    ##################
    #--- ANALYSIS ---#
    ##################
    call analysis{
        input :
        results = merge_saige.merged_results,
        dict = merge_matrix.gene_json,
        lof = lof,
        docker = docker,
        prefix = prefix
     }
}

task analysis {

    File results
    File dict		

    String? analysis_docker
    String docker
    String? final_docker = if defined(analysis_docker) then analysis_docker else docker    

    String lof
    String prefix
    Int disk_size = ceil(size(results,'GB')) + 10
    Int mem = ceil(size(results,'GB')) + 4
    
    String out_json =  "${prefix}_${lof}_gene_results.json"
    
    command <<<
     python3 /Scripts/analysis_results.py \
    -o . \
    --prefix "${prefix}_${lof}" \
    --saige-file ${results} \
    --dict-file ${dict}

    mv ${dict} ${out_json}
    >>>
    output {
        File json = out_json
        File saige_results = sub(out_json,'.json','.txt.gz')
        File pdf = sub(out_json,'.json','.pdf')
        }

     runtime {
        docker: "${final_docker}"
        cpu: 1
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b"
        preemptible: 1
        noAddress: true
    	}

}


task merge_saige {

    Array[File] gene_results
    
    String? merge_docker
    String docker
    String? final_docker = if defined(merge_docker) then merge_docker else docker    

    Int disk_size = ceil(size(gene_results[0],'GB')*length(gene_results)) + 10

    command <<<
    echo ${disk_size}
    # merge the results
    head -n1 ${gene_results[0]} > results.txt
    while read f; do echo $f && cat $f | sed -E 1d >> results.txt ; done < ${write_lines(gene_results)}
    gzip results.txt

    >>>

    output {
        File merged_results = "results.txt.gz"
        }

     runtime {
        docker: "${final_docker}"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b"
        preemptible: 1
        noAddress: true
    	}

    }

task pheno_saige {

    File gene_vcf
    File gene_tbi

    String pheno
    String variance_basepath 
    File variance = sub(variance_basepath,"PHENO",pheno) 
    File rda= sub(variance, ".varianceRatio.txt", ".rda")

    File samples
    String docker
    
    String outfile =  pheno +  ".SAIGE.txt"
    String disk_size = ceil(size(gene_vcf,'GB')) + 10
    
    command <<<
    step2_SPAtests.R \
    --vcfFile=${gene_vcf} \
    --vcfFileIndex=${gene_tbi} \
    --GMMATmodelFile=${rda} \
    --varianceRatioFile=${variance} \
    --sampleFile=${samples} \
    --SAIGEOutputFile=tmp.txt \
    --IsOutputAFinCaseCtrl=TRUE \
    --LOCO=FALSE \
    --IsOutputNinCaseCtrl=TRUE \
    --chrom 1

    #PREPEND PHENO TO COLUMNS AND REPLACE `SNPID` WITH `GENE`
    head -n1 tmp.txt | sed 's|SNPID|GENE|g' |  awk '{print "PHENO "$0}' > ${outfile}
    cat tmp.txt  | sed -E 1d |  awk '{print "${pheno} "$0}' >> ${outfile}
    >>>

    output {	
	File out = outfile
	}
	
    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    } 

}

task gene_matrix {

    String chrom_file
    String chrom
     
    File chrom_vcf = sub(chrom_file,'CHROM',chrom)
    File chrom_tbi = chrom_vcf + ".tbi"
    File lof_variants
    File samples  

    Float call_filter
    String call_type = if call_filter >0 then "--hard-call " + call_filter else "--gp0"
    
    Float? info_score
    String info_filter = if defined(info_score) then "--info_score " + info_score else " "
    Float max_maf
    
    String lof          
    String? matrix_docker
    String docker
    String? final_docker = if defined(matrix_docker) then matrix_docker else docker
    Int cpu
    String disk_size = 2*ceil(size(chrom_vcf,'GB'))+10
    
    command{
    echo ${disk_size}
    python3 /Scripts/lof.py \
        --lof_variants ${lof_variants} \
        --vcf ${chrom_vcf} \
        -o /cromwell_root/ \
        -c ${chrom} \
        ${info_filter} \
        --lof ${lof} \
        ${call_type} \
        -s ${samples} \
        --cpus ${cpu} \
        --max-maf ${max_maf} 
     }

     runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
         disks: "local-disk " + "${disk_size}" + " HDD"
        zones: "europe-west1-b"
        preemptible: 1
	bootDiskSizeGb: 20
	memory: "8 GB"
    }
    output {	
        File vcf = "/cromwell_root/${chrom}/${lof}_${chrom}_gene.vcf.gz"
        File gene_dict = "/cromwell_root/${chrom}/${lof}_${chrom}_gene_variants_dict.txt"
       }

}

task merge_matrix {

    Array[File] vcfs
    Array[File] dicts
    
    String prefix
    String lof

    String? merge_docker
    String docker
    String? final_docker = if defined(merge_docker) then merge_docker else docker

    Int disk_size = ceil(size(vcfs[0],'GB')*length(vcfs)) + 10

    String out_vcf = "${prefix}_${lof}_gene.vcf.gz"
    String out_json = "${prefix}_${lof}_gene.json"

    command <<<
    echo ${disk_size}
    # MERGE VCFS
    bcftools concat -f ${write_lines(vcfs)} -Oz -o ${out_vcf}
    tabix -f ${out_vcf}

    # MERGE JSON FILES FROM ALL CHROMS
    touch ${out_json} 
    while read f; do echo $f && jq -n 'reduce inputs as $i ({}; . * $i)' ${out_json} $f > tmp && mv tmp ${out_json}; done < ${write_lines(dicts)}
    
    >>>
     
     runtime {
        docker: "${final_docker}"
        cpu: "1"
	disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b"
        preemptible: 1
	bootDiskSizeGb: 20
	memory: "4 GB"
    }

    output{
        File gene_vcf = out_vcf
        File gene_tbi = "${out_vcf}.tbi"
        File gene_json = out_json
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
