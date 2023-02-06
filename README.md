# LOFanalysis

Regenie pipeline for LOF.
Burden mode for regenie is used. This pipeline does all the work from annotation/vcfs to results. 

# WDL
The wdl(wdl/regenie_lof.wdl) runs the whole pipeline.
The first step is the extraction of LOF variants. The task `extract_variants` does the job. It also filters out common variants (`max_maf` input) and low quality imputed variants (`regenie_lof.extract_variants.info_filter`)

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
