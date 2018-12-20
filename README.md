# LOFanalysis

This script generates the lof variants and outputs the gene_to_sample matrix based on the info score of the batch

## LOF.py

### Inputs:
`--annotated_file ` tsv.gz file with columns named gene and `$LOF`. It's used for mapping a variant to gene and LOF\
`--lof` : type of LOF. At the moment it accepts `most_severe` and `hc_lof`\
`-o` : out path\
`--bed` : path to plink bed file \
`--exclude `: File(s) with list of variants to exclude\
`--remove ` : File with list of samples to remove.\
`--samples ` : File with list of samples to use.\
`--cpus `: Number of parallel processes to run, by default the number of cpus of the machine\
`--maxMAF `: maximum minor allele frequency \
`--pargs` : the arguments to pass to the plink write-snplist call \
`--test`  : accepts an integer which is the number of genes to run per cpu. By deafult 0, which runs all genes 

E.g. `python3 ./Scripts/LOF.py --annotated_file /mnt/disks/tera/Data/R2_vep_annotated.tsv.gz --lof most_severe -o /mnt/disks/tera/LOF_test --bed /mnt/disks/tera/LOF/plink_test/most_severe.bed --exclude /mnt/disks/tera/Data/variants/lq_variants_0.9.txt /mnt/disks/tera/Data/variants/r2_blacklist_all.txt --samples /mnt/disks/tera/Data/R2_final_samples.txt --pargs "--maf 0.05"`

### How it works

The script first reads through the annotated file and saves the LOF carry variants, according to the type of LOF requested. It also builds a variant_to_gene dict.\
The final snplist is created with plink, passing the LOF variants, the variants to be removed and is then used to build a new small plink file, which reorders the samples to match the sample list passed.\
From there  the `.raw` matrix is created, which has the sample to variant LOF matrix.\
Using the list of snpslist and the variant_to_gene_dict a reverse gene_to_variants dict is built so that variants can be merged into genes. The genelist is split into `$cpus` chunks (a smaller subset if `$test` is used) to be run in a different cpu. Each process produces a temporary gene to sample matrix. The matrices are then pasted together adding a top row for finngen ids.\

## Docker

The `Docker` folder contains the script `build_docker.py` which builds the docker. `sudo python3 build_docker.py --version 0.001` is the command to run. 
Within my machine (pete-pet) one can then run `sudo docker run -v /mnt/disks/tera/:/mnt/disks/tera/ -it eu.gcr.io/finngen-refinery-dev/lof:0.001 /bin/bash` to test it locally in order to avoid having to download the massive R2 plink files.
Once the docker is tested, it's enough to rerun `build_docerk.py` with the `--push` flag and it will be pushed to gcloud.

## wdl

`saige.wdl` is the wdl to run. Nothing special to note except for a couple of hacks:
- The full/complete version of the wdl reuqires the R2 plink file as an input. For testing purposes in case of a new annotated file, it's possible to cut the times by using the plink files generated locally by the `LOF.py` script. In this way the matrix calculation will be extremely fast (10 minutes or so) and the output will be identical. 
- The flag `--pargs` of the python script is a hack to allow to pass more flags to plink withouth having to modify the python script specifically. However, the pargs are passed only to the part of code which outputs the snplist from the original plink file, so it does not affect the creation of the matrix.

