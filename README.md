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
FG2A3ANP94 0.00199711 0.00248884 0.00362434 -0.00017248 0.00528694 0.00986158 0.00511827 -0.00268237 0.00628966 0.00071082
FG2A3QFM2N 0.000884429 0.00301878 -0.00271975 -0.00184375 0.00108309 0.00694645 -0.00106438 0.00211031 -0.00620706 0.00118039
FG2A3WKJV9 -0.00287612 0.000871676 0.00186428 -0.000480518 -0.00661943 0.00374202 -0.00737054 -0.00769109 -0.00201046 0.00240494
FG2A4XB9BM 0.0105941 -0.00832174 0.00698353 -0.00713543 -0.000172974 -0.000327632 0.00490943 0.0011419 -0.00601704 -0.0112859
FG2A5JWZXA 0.0026648 -0.00379872 -0.00824472 -0.00314098 0.00302931 -0.00361574 -0.0207097 0.0074027 -0.0078156 -0.00118258
FG2A6HHSVG -0.000965202 0.00869009 -4.33087e-05 0.00257666 0.00617113 0.000641854 0.00157645 0.00142078 0.00672996 -0.00481868
FG2A7RYAVJ 0.00599325 0.00179053 0.00319636 0.00105346 0.000741381 -0.00167067 -0.00577088 -0.00310104 0.000566559 0.00197225
FG2A7X9A7Q 0.00578686 -0.00145477 0.00194329 -0.0124015 -0.00525735 -0.00217946 -0.00206569 0.00202294 0.00300609 0.000458245`
* a phenotype file, where the first columns are samples and the columns are the various phenotypes 
* pheno-list.txt a file with the list of phenotypes that need to be run 