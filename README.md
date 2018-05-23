# LOFanalysis

## filter_variants.py

This script generates the lof variants and outputs the gene_to_sample matrix based on the info score of the batch

### filter the variants
The first step is to generate the variants based on the annotated File.

The comand is \
`python filter_variants.py filter --annotatedFile ../Data/annotated_variants.gz --lof hc_lof`

where `annotateFile` is a file where a column is marked as `lof`.

### generate the matrix
Once the variants are generated one needs to run the wdl script to filter out the variants from the original data. Once the data is generated one needs to run the `generate-matrix` command\
`python filter_variants.py generate-matrix --lof hc_lof --plinkPath ~/hc_lof_data/ --oPath ~/results/`

where\
`plinkPath` is the path to the plink binaries
`oPath` is the path where the results are ouputted to. It will automatically be changed to oPath + lof.

The result is a file called `$lof_gene_to_sample.tsv` matrix where each row is the lof data of a gene.

## gene_analysis.py

This script producedes the logit regression on data based on phenotypes.

For this run we need the following data:\
* the pc data, organized in a file where the first columns are samples and the rest are space separated floats.
* a phenotype file, where the first columns are samples and the columns are the various phenotypes 
* pheno-list.txt a file with the list of phenotypes that need to be run 