# LOFanalysis

This script generates the lof variants and outputs the gene_to_sample matrix based on the info score of the batch

## LOF.py
Usage:
```
usage: lof_test.py [-h]
                   (--annotation ANNOTATION | --lof_variants LOF_VARIANTS) -o
                   OUT_PATH -c CHROM --vcf VCF --lof {hc_lof,most_severe}
                   [--info_score INFO_SCORE] [-s SAMPLE_FILE] [--cpus CPUS]
                   [--force] [--test]
```
### Inputs:
Required:\
One between:\
`--annotation ` tsv.gz file with columns named gene and `$LOF`. It's used for mapping a variant to gene and LOF\
`--lof_variants` tsv file where the first column is the variant and the second is the gene. It's the output of the previous flag. \

`--lof` : type of LOF. At the moment it accepts `most_severe` and `hc_lof`\
`-o` : out path\
`--vcf` : path to vcf file \
`-c` : chromsome number 


Optional:\
`--samples ` : File with list of samples to use.\
`--cpus `: Number of parallel processes to run, by default the number of cpus of the machine\
`--test`  : accepts an integer which is the number of genes to run per cpu. By deafult 0, which runs all genes 


### How it works

The script first reads through the annotated file and saves the LOF carry variants, according to the type of LOF requested. It also builds a variant_to_gene dict.\
Then the variants are split into chunks for speedup purposes. For each variant chunk a `bcftools` command filter the vcf for those positions and returns the probability of GP=0 for the variant in a variant to sample matrix.\
Finally, the variant chunks are merged and transposed so a final sample to variant to matrix is built.\
Now variants are merged into genes with the sample to gene value being 1 minus the product of the GP=0 for each variant to obtain the final gene to sample variant.

## Docker

The `./Docker` folder contains the script `build_docker.py` which builds the docker. `sudo python3 build_docker.py --version $VERSION` is the command to run it. 
Once the docker is tested, it's enough to rerun `build_docker.py` with the `--push` flag and it will be pushed to gcloud.

## wdl


The wdl runs the following steps:
1) it generates the LOF matrix for each chromosomes 
2) it merges the matrix into one single files 
3) it extracts the phenotype names from the null data 
4) it runs SAIGE for each pheno 
5) it merges the results and sorts them by pval