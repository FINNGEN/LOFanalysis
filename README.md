# LOFanalysis

## filter_variants.py

This script generates the lof variants and outputs the gene_to_sample matrix based on the info score of the batch

The first step is to generate the variants based on the annotated File.

The comand is \
`python filter_variants.py filter --annotatedFile ../Data/annotated_variants.gz --lof hc_lof`

where `annotateFile` is a file where a column is marked as `lof`.