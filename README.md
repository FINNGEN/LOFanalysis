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

`lof_list` determines which annotations are incluced for flagging variants as LOF.
`info_filter` puts a minimimum threshold of info score across all batches.
`max_max` filters out common variants above such threshold


Once we have a list of LOF variants, the wdl proceeds to subet the input vcfs, creating a merged vcf as well as a merged bgen.