workflow filter_lof {

	Array[String] chromList
	
	scatter (chrom in chromList) {
		call chrom_filter_merge {
		     input:
			chrom = chrom
	     }
	}

	call fix_files {
	     input:
		all_files = chrom_filter_merge.all_files
	}

	output {	       
	       Array[Array[File]] output_files = chrom_filter_merge.all_files
	       Array[Float] freedisk = chrom_filter_merge.free_disk
	}

}

task chrom_filter_merge {

     String chrom
     String chromPath
     String tbiPath
     String variantsPath
     
     File cFile = sub(chromPath,"_CHROM",chrom)
     File tbiFile = sub(tbiPath,"_CHROM",chrom)
     File vFile = variantsPath	
     
     

     command {    	     
          python3 /ProcessImputedData/Scripts/filter_tabix.py --cFile ${cFile} --tbiFile ${tbiFile} --vFile ${vFile} --chrom ${chrom} --oPath "/cromwell_root/results/" --separator _
	  du -sh | cut -f1 | rev | cut -c 2- | rev > freedisk.txt
	  } 

     runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/filter_variants:0.002"
	cpu: 3
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b"
        preemptible: 0
	bootDiskSizeGb: 20

    }

   


   output {
       Array[File] all_files = glob("results/chrom_merge/chr*") 
       Float free_disk = read_float("freedisk.txt")
   }

}
   

task fix_files {

     Array[Array[File]] all_files
     
     String output_file = "chrom_merged"
     
     command {
     	     python3 /ProcessImputedData/Scripts/fix_files.py --oFile ${output_file} --oPath "/cromwell_root/results/" ${sep=' ' all_files}
     }

     runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/filter_variants:0.002"
	cpu: 1
        disks: "local-disk 1000 HDD"
        zones: "europe-west1-b"
        preemptible: 0
	bootDiskSizeGb: 20
	memory: "4 GB"

    }

    output {
        Array[File] ouput_files = glob("results/" + output_file + "*")
	}
}



