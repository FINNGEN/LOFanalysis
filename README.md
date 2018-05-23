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
* the pc data, organized in a file where the first columns are samples and the rest are space separated floats\eg.
`FG2A2F47CU 0.00559777 0.000546521 0.00166026 -0.00227628 -0.00182399 0.00103257 0.000670258 -0.00127393 -0.00170497 -6.6967e-05
FG2A2V36PK -0.00338977 -0.00834218 -0.00990391 -0.00138755 -0.00604381 1.23167e-05 0.00367658 -0.00100084 -0.00268119 0.00170265
FG2A3ANP94 0.00199711 0.00248884 0.00362434 -0.00017248 0.00528694 0.00986158 0.00511827 -0.00268237 0.00628966 0.00071082`
* a phenotype file, where the first columns are samples and the columns are the various phenotypes 
* pheno-list.txt a file with the list of phenotypes that need to be run 