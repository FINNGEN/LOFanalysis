# LOFanalysis

Regenie pipeline for LOF.
Burden mode for regenie is used. This pipeline does all the work from annotation/vcfs to results. 

# WDL
The wdl(wdl/regenie_lof.wdl) runs the whole pipeline.
The first step is the extraction of LOF variants. The task `extract_variants` does the job. It also filters out common variants (`regenie_lof.max_maf` input) and low quality imputed variants (`regenie_lof.extract_variants.info_filter`)

```
chr10_1020846_AG_A	IDI2	frameshift_variant
chr10_102733581_C_T	SFXN2	stop_gained
chr10_103090970_C_CCAGA	NT5C2	frameshift_variant
chr10_103246174_C_T	RPEL1	stop_gained
chr10_103419670_C_T	PDCD11	stop_gained
chr10_103449599_CAG_C	CALHM2	frameshift_variant
chr10_103455672_G_A	CALHM1	stop_gained
chr10_103473333_C_T	CALHM3	stop_gained
chr10_113595655_G_A	NRAP	stop_gained
chr10_113597145_TG_T	NRAP	frameshift_variant
```

The variants are then passed to `convert_vcf` and then to `merge` to build a single bgen with only lof varirants.

From there a scatter takes each pheno and runs the `regenie` task, where the association takes place.

The task `merge_results` then produces all files needed for release:
- a readme file
- all hits sorted by logp
- merged logs
- a summary of the results only including sig hits and sorted by difference between mlgop and max(mlogp) of the variants in the gene
- the file required for sql import

##

Relevant inputs:

`extract_variants.annot_file`: File with lof, AF and info score data. All information is extracted automatically from header names.'\
`merge.bargs`: the vcf--> bgen conversion params \
`extract_variants.lof_list`: the list of lof keywords to include \
`regenie.cov_file`: self explanatory, where the covariate and pheno data is located \
`regenie.masks_type`: the bgen parameter that determines the merging logic for variants into genes \
`regenie.bargs` : all other regenie args.