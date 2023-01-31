version development

workflow regenie_lof {

  input {
    String docker
    String prefix
    Float max_maf
    Array[String] chrom_list
  }

  # get lof variants
  call extract_variants { input: docker = docker, max_maf = max_maf}
 
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

  zcat ~{annot_file} | awk -v OFS='\t' -v c1=$AIND -v c2=$IIND -v c3=$GIND -v c4=$MIND '{print $1,$c1,$c2,$c3,$c4}' | awk '$2 < ~{max_maf} || $2 > ~{1-max_maf}' | awk '$3 > ~{info_filter}' |  grep -wf ~{write_lines(lof_list)} |  tr ':' '_' | awk '{print "chr"$0}' | cut -f 1,4,5 > lof_variants.txt
  
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
    }

}
