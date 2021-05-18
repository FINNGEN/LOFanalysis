import argparse,gzip,os
import numpy as np
from file_utils import file_exists,make_sure_path_exists
from collections import defaultdict


def annotation_filter(annot_file,variant_file,out_annot):
    """
    Function that filters the small annotation file to only include variants in the lof pgen.
    """
    variants = set(np.loadtxt(variant_file,dtype = str,usecols = 0))
    with gzip.open(annot_file) as i,open(out_annot,'wt') as o:
        next(i)
        for line in i:
            variant,_,gene,most_severe = line.strip().decode("utf-8").split()
            variant = "chr" + "_".join(variant.split(':'))
            if variant in variants:
                o.write("\t".join([variant,gene,most_severe]) + "\n")


def create_sets(annot,set_file):
    """
    Builds sets for regenie

    """
    gene_dict = defaultdict(list)

    with open(annot) as i:
        for line in i:
            snp,gene,most_severe = line.strip().split()
            gene_dict[gene].append(snp)

    print('gene dict built')
    with open(set_file,'wt') as o:
        for gene in gene_dict:
            variants = gene_dict[gene]
            chrom_list,pos_list = zip(*[elem.split("_")[:2] for elem in variants])
            chrom_set = set(chrom_list)
            if len(chrom_set) ==1:
                chrom = chrom_list[0].split("chr")[1]
                o.write('\t'.join([gene,chrom,pos_list[0],','.join(variants)]) +'\n')
            else:
                print(gene,chrom_set)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Regenie LOF annotation")
    parser.add_argument('-o',"--out",type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--annot",type = file_exists,help = "annotation file",required = True)
    parser.add_argument("--variant-file",type = file_exists,help = "annotation file",required = True)

    args=parser.parse_args()
    make_sure_path_exists(args.out)

    out_annot = os.path.join(args.out,"annot.tsv")
    out_set = os.path.join(args.out,"sets.tsv")

    annotation_filter(args.annot,args.variant_file,out_annot)
    create_sets(out_annot,out_set)
