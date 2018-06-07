# LOFanalysis

## filter_variants.py

This script generates the lof variants and outputs the gene_to_sample matrix based on the info score of the batch

### filter the variants
The first step is to generate the variants based on the annotated File.

The comand is \
`python filter_variants.py --annotatedFile ../Data/annotated_variants.gz --lof hc_lof filter`

where `annotatedFile` is a file where a column is marked as `$lof`. In case of `most_severe`, it searches for a list of other labels in the lof column.

### generate the matrix
Once the variants are generated one needs to run the wdl script to filter out the variants from the original data. Once the data is generated one needs to run the `generate-matrix` command\
`python filter_variants.py --lof hc_lof --annotatedFile ../Data/annotated_variants.gz generate-matrix --plinkPath ~/Data/hc_lof_data/chrom_merged --oPath ~/results/ --samplePath ~/LOFanalysis/Data/sample_info.txt`

where:\
`plinkPath` is the path to the plink binaries\
`oPath` is the path where the results are ouputted to. It will automatically be changed to oPath + lof.\
`samplePath` is the path to the sample to batch mapping. It needs to be formatted SAMPLE:BATCH, with possible other column separated elements, as long as sample is the first one and batch the last.\

The result is a file called `$lof_gene_to_sample.tsv` matrix where each row is the lof data of a gene. The values are the maximum info score value across all variants that belong to the gene, where the info score is the batch info score for that variant and sample.

## gene_analysis.py

At the moment this script is not implemented.