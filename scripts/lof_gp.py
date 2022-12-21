#simplified version where the input is already the vcf filtered with only LOF. SO the input should be only the mapping really and it should be processed as a whole. I should try to loop over the genes directly, removing a lot of intermediate steps.

from collections import defaultdict
from itertools import chain
from functools import partial 
import numpy as np
import argparse,os,multiprocessing,logging,time
from file_utils import file_exists,make_sure_path_exists,tmp_bash,spin_bash,progressBar,pretty_print,log_levels,mapcount,cpus
import pandas as pd

######################################
#------------VARIANTS----------------#
######################################
def gene_lof_map(lof_map):
    """
    Returns gene to variant dict sorted by count.
    """
    lof_dict = defaultdict(list)
    with open(lof_map) as i:
        for line in i :
            variant,gene,*_ = line.strip().split()
            lof_dict[gene].append(variant)

    return lof_dict

def extract_vcf_variants(vcf,lof_dict,out_path):
    """
    Check that all variants are in mapping
    """
    vcf_tmp = os.path.join(out_path,"vcf_variants.txt")
    if not os.path.isfile(vcf_tmp):
        print('Extracting variants from vcf...')
        cmd = f"zcat {vcf} | grep -v '#' |  cut -f 3 > {vcf_tmp}"
        spin_bash(cmd)
    return np.loadtxt(vcf_tmp,dtype=str)

def shared_variants(vcf,lof_map,out_path):

    lof_dict = gene_lof_map(lof_map)
    logging.info(f"{len(lof_dict.keys())} genes in the mapping.")
    lof_variants = list(chain(*lof_dict.values()))
    logging.info(f"{len(lof_variants)} variants across all genes.")

    vcf_variants = extract_vcf_variants(vcf,lof_dict,out_path)
    logging.info(f"{len(vcf_variants)} variants in the vcf file.")

    # filter out missing variants in the vcf.
    new_dict = defaultdict(list)
    for gene in lof_dict:
        variants =[elem for elem in lof_dict[gene] if elem in vcf_variants]
        if variants:
            new_dict[gene] = variants
            
    ordered_genes = sorted(new_dict, key=lambda k: len(lof_dict[k]), reverse=True)
    logging.info(f"{len(ordered_genes)} genes left.")
    lof_variants = list(chain(*new_dict.values()))
    logging.info(f"{len(lof_variants)} variants across all genes.")
    with open(os.path.join(out_path,'gene_dict.tsv'),'wt') as o:
        for gene in new_dict:
            o.write(gene + '\t' + ','.join(new_dict[gene]) + '\n')
            
    return new_dict,ordered_genes


#####################################
#------------  VCF  ----------------#
#####################################

def merge_gene_chunks(lof_genes,lof_dict,vcf,out_path,test):
    """
    Creates the vcf from the gene chunks
    """

    sample_header= os.path.join(out_path,'samples_header.txt')
    if not os.path.isfile(sample_header):
        logging.info("Extracting samples")
        sample_cmd = f" echo VARIANT > {sample_header} && bcftools query -l {args.vcf} test_filter >> {sample_header}"
        spin_bash(sample_cmd)

    tmp_path = os.path.join(out_path,'tmp')
    make_sure_path_exists(tmp_path)
    
    # MULTIPROC APPROACH
    genes = lof_genes if not test else lof_genes[:cpus]
    print(f"{len(genes)} genes to parse.")
    print(f"{cpus} cpus being used.")
    pool = multiprocessing.Pool(cpus)
    multiproc_func = partial(gene_chunk,lof_dict=lof_dict,vcf=vcf,out_path=tmp_path,sample_header=sample_header)
    map(multiproc_func,genes)
    results = pool.map_async(multiproc_func,genes,chunksize=1)
    while not results.ready():
        time.sleep(2)
        progressBar(len(genes) - results._number_left,len(genes))
    progressBar(len(genes),len(genes))
    print('\ndone')
    #multiproc_func(lof_genes[10])
    
def gene_chunk(gene,lof_dict,vcf,out_path,sample_header):
    """
    Processes the gene as a single chunk. Then the chunk is transposed ( so that i can read in variants as columns)
    """
    logging.debug(gene)
    gene_chunk = os.path.join(out_path,gene + '.tmp')
    logging.debug(gene_chunk)
    chrom = lof_dict[gene][0].split('_')[0]
    #use bcftools to keep positions
    gene_variants = os.path.join(out_path,gene + '.variants')
    with open(gene_variants,'wt') as o:
        for variant in lof_dict[gene]:o.write(variant + '\n')
    # WRITE SAMPLES HEADER
    header_cmd = f"cat {sample_header}| tr '\\n' '\\t' > {gene_chunk} "
    # TREAT VARIANTS <0.5
    variant_filter = f" -i 'ID=@{gene_variants} & AF <= 0.5 ' "
    min_af_cmd = f"bcftools query --regions {chrom}  -f " + "'%ID[\\t%GP{0}]\\n' " + f" {variant_filter} {vcf} >> {gene_chunk} "
    # TREAT VARIANTS >0.5
    variant_filter = f" -i 'ID=@{gene_variants} & AF > 0.5 ' "
    max_af_cmd = f"bcftools query --regions {chrom}  -f " + "'%ID[\\t%GP{2}]\\n' " + f" {variant_filter} {vcf} >> {gene_chunk} "
    # JOIN COMMANDS
    matrix_cmd = f"{header_cmd} && {min_af_cmd} && {max_af_cmd}"
    logging.debug(matrix_cmd)
    tmp_bash(matrix_cmd)

    logging.debug('Reading in data and merging')
    merged_gene = os.path.join(out_path,gene + '.txt')
    logging.debug(merged_gene)
    data = 1 - pd.read_csv(gene_chunk,sep ='\t',index_col=0).T.prod(axis = 1).values
    with open(merged_gene,'wt') as o : o.write(gene +'\t'.join(list(map(str,data)))+'\n')
        
def main(args):

    pretty_print("PROCESSING VARIANTS")
    lof_dict,ordered_genes = shared_variants(args.vcf,args.lof_map,args.out_path)

    pretty_print("GENES")
    merge_gene_chunks(ordered_genes,lof_dict,args.vcf,args.out_path,args.test)
                             
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Builds vcf with custom processing of gene merging.")
    parser.add_argument('-o',"--out_path",type = str, help = "folder in which to save the results", required = True)
    # MAPPING
    parser.add_argument("--lof_map", type= file_exists,help="Chrom to gene mapping for lof variants")
    # CHUNK CREATION
    parser.add_argument("--vcf", type= file_exists,help="path to vcf file",required = True)
    # MERGING
    #call_parser.add_argument("--gp0", action='store_true', help =  "Works with GP = 0")
    #MISC
    parser.add_argument( "-log",  "--log",  default="warning", choices = log_levels, help=(  "Provide logging level. " "Example --log debug', default='warning'"))
    parser.add_argument('--test',action = 'store_true',help = 'Runs test version')

    args=parser.parse_args()
    # logging level
    level = log_levels[args.log]
    logging.basicConfig(level=level,format="%(levelname)s: %(message)s")

    make_sure_path_exists(args.out_path)
    main(args)
