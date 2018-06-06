#!/bin/bash


python filter_variants.py --lof hc_lof generate-matrix --plinkPath ~/Data/hc_lof_data/chrom_merged --oPath ~/results/ --samplePath ~/LOFanalysis/Data/sample_info.txt 
python filter_variants.py --lof most_severe generate-matrix --plinkPath ~/Data/most_severe_data/chrom_merged --oPath ~/results/ --samplePath ~/LOFanalysis/Data/sample_info.txt 
