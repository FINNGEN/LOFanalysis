# LOFanalysis

Regenie pipeline for LOF.

# BASIC DATA MUNGING 


The [filter_lof wdl](./wdl/filter_lof.wdl) filters the input vcfs to only use LOF variants. The input variants are based on file annotations using VEP.

These are the parameters in the `json` file that are relevant to the the filtering.

```
"filter_lof.extract_variants.lof_list":  ["frameshift_variant","splice_donor_variant","stop_gained","splice_acceptor_variant"],
"filter_lof.extract_variants.info_filter": "0.95",
"filter_lof.extract_variants.max_maf": "0.05",
```

`lof_list` determines which annotations are incluced for flagging variants as LOF. \
`info_filter` puts a minimimum threshold of info score across all batches. \
`max_max` filters out common variants above such threshold.


Once we have a list of LOF variants, the wdl proceeds to subet the input vcfs, creating a merged vcf as well as a merged bgen and a summmary file with the variant to gene mapping of containing all the variants used.


## GP
This is the explanation of the GP based LOF analysis.

The [gp_lof wdl](./wdl/gp_lof.wdl) in this case does the heavy lifting and is based on the [lof_gp.py](./scripts/lof_gp.py) python script. The script takes as an input a variant to gene tsv file like the one produced by the `filter_lof.wdl`
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
and for each gene in the set it takes all variants in the input vcf and merges them to create a custom dosage value. The logic being is that the dosage of the gene represents the probability of carrying at least one LOF variants based on the genotype GP for each variant.
For sample $i$ and gene $g$ the custom dosage would be:\
$$D_{g}^{i} = 1 - \prod_{v \in g} GP_{v}^{i}[maj]$$ where $GP[maj]$ depends on the AF of the variant.