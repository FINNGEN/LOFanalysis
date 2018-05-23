# LOFanalysis

## filter_variants.py

This script generates the lof variants and outputs the gene_to_sample matrix based on the info score of the batch

The first step is to generate the variants based on the annotated File.

The comand is \
`python filter_variants.py filter --annotatedFile ../Data/annotated_variants.gz --lof hc_lof`

where `annotateFile` is a file where a column is marked as `lof`.

Once the variants are generated one needs to run the wdl script to filter out the variants from the original data. Once the data is generated one needs to run the `generate-matrix` command\
`python filter_variants.py generate-matrix --lof hc_lof --plinkPath ~/hc_lof_data/ --oPath ~/results/`

where\
`plinkPath` is the path to the plink binaries
`oPath` is the path where the results are ouputted to. It will automatically be changed to oPath + lof.

## gene_analysis.py

