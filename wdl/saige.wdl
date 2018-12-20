task pheno_saige {
	
	File nullfile
	File lofMatrix
	File varianceratiofile = sub(nullfile, ".rda", ".varianceRatio.txt")
	File samplefile
	Int minmac
	String docker
	Int cpu
	Float mem
	String loco
	String outfile = basename(nullfile, ".rda") +  ".SAIGE.txt"

	command {

	step2_SPAtests.R \
	--dosageFile=${lofMatrix} \
	--dosageFileNrowSkip=1 \
	--dosageFileNcolSkip=1 \
	--dosageFilecolnamesSkip="GENE" \
	--minMAC=10 \
	--sampleFile=${samplefile} \
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


workflow LOF_saige{

	File null_list

	Array[String] nullfiles = read_lines(null_list)

	scatter (nullfile in nullfiles){
	    call pheno_saige{
	    	 input : nullfile=nullfile
		 }
	}
}