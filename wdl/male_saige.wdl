task pheno_saige {
	
	
	File matrix
	File varianceratiofile
	File nullfile= sub(varianceratiofile, ".varianceRatio.txt", ".rda")
	File samples
	Int minmac
	String docker
	Int cpu
	Float mem
	String loco
	String outfile = basename(nullfile, ".rda") +  ".MALE.SAIGE.txt"

	command {

	step2_SPAtests.R \
	--dosageFile=${matrix} \
	--dosageFileNrowSkip=1 \
	--dosageFileNcolSkip=1 \
	--dosageFilecolnamesSkip="GENE" \
	--minMAC=10 \
	--sampleFile=${samples} \
	--LOCO=${loco} \
	--numLinesOutput=1000 \
	--IsOutputAFinCaseCtrl=TRUE \
	--GMMATmodelFile=${nullfile} \
	--varianceRatioFile=${varianceratiofile} \
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




task lof_matrix {


     File bed_file	
     File fam_file = sub(bed_file,".bed",".fam")
     File bim_file = sub(bed_file,".bed",".bim")
     
     File annotated_file
     File exclusion_file
     Array[File] exclusion_files = read_lines(exclusion_file)
     File samples
     
     String? pargs
     String plink_args = if defined(pargs) then " --pargs " + pargs else ""
     String LOF
     
     String docker
     String cpu
     String mem
     String disk_size
     String out  ="/cromwell_root/results/"
     
     command {

     python3 /Scripts/LOF.py \
     --annotated_file ${annotated_file} \
     --lof ${LOF} \
     -o ${out} \
     --bed ${bed_file} \
     --exclude ${sep=' ' exclusion_files }  \
     --samples ${samples} \
     ${plink_args}    

     }

     output {
     File matrix = "${out}" + "${LOF}" + "_matrix.txt"
     File gene_dict = "${out}" + "${LOF}" + "_gene_variants_dict.txt"
  }


     runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b"
        preemptible: 2
	
    }

     }


workflow LOF_saige{

 
	File male_samples
		
	call lof_matrix {
	     input: samples = male_samples
		        }

	Int minmac
	String docker
	Int cpu
	Float mem
	String loco


	File male_variance_list
	Array[String] male_variances = read_lines(male_variance_list)
	
	scatter (male_variance in male_variances){
	    call pheno_saige{
	    	 input :
		       samples = male_samples,
		       varianceratiofile= male_variance,
		       matrix = lof_matrix.matrix,
		       minmac = minmac,
		       docker = docker,
		       cpu= cpu,
		       mem=mem,
		       loco =loco
		 }
	}


}