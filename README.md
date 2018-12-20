# LOFanalysis

This script generates the lof variants and outputs the gene_to_sample matrix based on the info score of the batch

## Inputs:
`--annotated_file ` tsv.gz file with columns named gene and `$LOF`. It's used for mapping a variant to gene and LOF\
`--lof` : type of LOF. At the moment it accepts `most_severe` and `hc_lof`\
`-o` : out path\
`--plink` : path to plink files (e.g. `/foo/bar/data`) where `data` is the root of the bim/bed/fam files\
`--exclude `: File(s) with list of variants to exclude\
`--remove ` : File with list of samples to remove`.\
`--samples ` : File with list of samples to use`.\
`--cpus `: Number of parallel processes to run, by default the number of cpus of the machine\
`--maxMAF `: maximum minor allele frequency \
`--pargs` : the arguments to pass to the plink write-snplist call \
`--test`  : accepts an integer which is the number of genes to run per cpu. By deafult 0, which runs all genes \

E,g, `python3 ./Scripts/LOF.py --annotated_file /mnt/disks/tera/Data/R2_vep_annotated.tsv.gz --lof most_severe -o /mnt/disks/tera/LOF_test -p /mnt/disks/tera/LOF/plink_test/most_severe --exclude /mnt/disks/tera/Data/variants/lq_variants_0.9.txt /mnt/disks/tera/Data/variants/r2_blacklist_all.txt --samples /mnt/disks/tera/Data/R2_final_samples.txt --pargs "--maf 0.05" `
`
## How it works

The script first reads through the annotated file and saves the LOF carry variants, according to the type of LOF requested. It also builds a variant_to_gene dict.\
Then, the final snplist is created with plink, passing the LOF variants, the variants to be removed and the samples that also need to be removed.\
Using the snplist, then a `.raw` matrix is created, which has the sample to variant LOF matrix.\
Also, the matrix will be rearranged so that the order is based on the list of samples provided, in case it's not alreaydy.\
Using the list of snpslist and the variant_to_gene_dict a reverse gene_to_variants dict is built so that variants can be merged into genes. The genelist is split into `$cpus` chunks (a smaller subset if `$test` is used) to be run in a different cpu. Each process produces a temporary gene to sample matrix. The matrices are then pasted together adding a top row for finngen ids.\




