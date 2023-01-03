#simplified version where the input is already the vcf filtered with only LOF. SO the input should be only the mapping really and it should be processed as a whole. I should try to loop over the genes directly, removing a lot of intermediate steps.

from collections import defaultdict
from itertools import chain
from functools import partial 
import numpy as np
import argparse,os,multiprocessing,logging,time,datetime
from file_utils import file_exists,make_sure_path_exists,tmp_bash,spin_bash,progressBar,pretty_print,log_levels,mapcount,cpus,int_or_float
import pandas as pd
from pathlib import Path

######################################
#------------VARIANTS----------------#
######################################
def gene_lof_map(lof_map):
    """
    Returns gene to variant dict sorted by count.
    """
    lof_dict = defaultdict(list)
    gene_chrom_map = {}
    with open(lof_map) as i:
        for line in i :
            variant,gene,*_ = line.strip().split()
            lof_dict[gene].append(variant)
            chrom = variant.split('_')[0].split('chr')[1]
            gene_chrom_map[gene] = int(chrom)
    
    return lof_dict,gene_chrom_map

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

    lof_dict,gene_chrom_map = gene_lof_map(lof_map)
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

    # order by chromosome
    ordered_genes = sorted(new_dict, key=lambda k: int(gene_chrom_map[k]))
    logging.info(f"{len(ordered_genes)} genes left.")
    lof_variants = list(chain(*new_dict.values()))
    logging.info(f"{len(lof_variants)} variants across all genes.")
    with open(os.path.join(out_path,'lof_gene_gp_dict.tsv'),'wt') as o:
        for gene in new_dict:
            o.write(gene + '\t' + ','.join(new_dict[gene]) + '\n')
            
    return new_dict,ordered_genes,gene_chrom_map


#####################################
#------------  VCF  ----------------#
#####################################

def merge_gene_chunks(lof_genes,lof_dict,vcf,out_path,test,force):
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
    gene_matrix = os.path.join(out_path,'lof_genes.txt')
    if not os.path.isfile(gene_matrix) or args.force:
        genes = lof_genes if not test else lof_genes[:cpus]
        print(f"{len(genes)} genes to parse.")
        print(f"{cpus} cpus being used.")
        pool = multiprocessing.Pool(cpus)
        multiproc_func = partial(gene_chunk,lof_dict=lof_dict,vcf=vcf,out_path=tmp_path,sample_header=sample_header)
        results = pool.map_async(multiproc_func,genes,chunksize=1)
        while not results.ready():
            time.sleep(2)
            progressBar(len(genes) - results._number_left,len(genes))
        progressBar(len(genes),len(genes))
        pool.close()
        print('\ndone')
        
        # dump results
        with open(gene_matrix,'wt') as o:
            for entry in results.get():
                o.write(entry + '\n')
    else:
        logging.info(f"{gene_matrix} already generated: using cache.")
        
    return gene_matrix,sample_header
    
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
    with open(gene_chunk,'wt') as o:o.write('\t'.join(np.loadtxt(sample_header,dtype=str)) + '\n')

    ## I need take into account allele frequency in order to choose which GP to extract.
    
    # TREAT VARIANTS <0.5
    variant_filter = f" -i 'ID=@{gene_variants} & AF <= 0.5 ' "
    min_af_cmd = f"bcftools query --regions {chrom}  -f " + "'%ID[\\t%GP{0}]\\n' " + f" {variant_filter} {vcf} >> {gene_chunk} "
    # TREAT VARIANTS >0.5
    variant_filter = f" -i 'ID=@{gene_variants} & AF > 0.5 ' "
    max_af_cmd = f"bcftools query --regions {chrom}  -f " + "'%ID[\\t%GP{2}]\\n' " + f" {variant_filter} {vcf} >> {gene_chunk} "
    # JOIN COMMANDS
    matrix_cmd = f"{min_af_cmd} && {max_af_cmd}"
    logging.debug(matrix_cmd)
    tmp_bash(matrix_cmd)

    logging.debug('Reading in data and merging')
    merged_gene = os.path.join(out_path,gene + '.txt')
    logging.debug(merged_gene)

    # This is where the action happens. The gene data is read merged.
    variant_data = pd.read_csv(gene_chunk,sep ='\t',index_col=0).T
    # data is aggregated across all variants to give us the probablity of carrying all the non LOF variants
    gene_data = variant_data.prod(axis=1).values
    # 1- gene_data gives us the probabily of carrying *at least* one LOF variant
    lof_data = np.round(1- gene_data,2)
    data = '\t'.join([gene] + list(map(str, lof_data)))
    
    return data


def build_vcf(gene_matrix,sample_header,gene_chrom_map,out_path):
    """
    Main function that builds the vcf. It reads in the chunk gene data and appends it to the vcf adding bogus pos/ref/alt data. Chrom instead is taken from the actual variant value.
    """
    # get vcf header
    parent_path = Path(os.path.realpath(__file__)).parent.parent
    vcf_header = os.path.join(parent_path,'data/','vcf_header.txt')

    # read in sample data
    samples = '\t'.join(np.loadtxt(sample_header,dtype=str)[1:])

    out_vcf = os.path.join(out_path,'lof_gene_gp.vcf')
    with open(out_vcf,'wt') as o,open(vcf_header,'rt') as i,open(gene_matrix) as gm:
        # update header info
        for line in i:
            line = line.replace('[DATE]',datetime.datetime.now().strftime("%Y%m%d"))
            if '[SAMPLES]' in line:
                line = line.replace('[SAMPLES]',samples)
            o.write(line)
            
        # loop through matrix and enter sample data
        ref=['','POS','ID','A','C','.','.','.','GT:GP']
        for j,line in enumerate(gm):
            # get gene data
            gene,*data = line.strip().split()
            chrom = gene_chrom_map[gene]
            logging.debug(f"{gene} {chrom}")
            # create bogus position and update data
            ref[0],ref[1],ref[2] = str(chrom),str(chrom*100000 + j),gene
            # go through entries and create bogus GTs and GPs
            gps = ["0|0:"+ ','.join(map(str,map(int_or_float,[round(1-float(elem),2),round(float(elem),2),0]))) for elem in data]
            line = '\t'.join(ref + gps) + '\n'
            o.write(line)

    logging.info("bgzipping...")
    cmd = f'bgzip -f  {out_vcf} '
    tmp_bash(cmd)
    logging.info("indexing ...")
    cmd = "tabix -f {out_vcf}.gz"
    tmp_bash(cmd)
    
def main(args):

    pretty_print("PROCESSING VARIANTS")
    lof_dict,ordered_genes,gene_chrom_map = shared_variants(args.vcf,args.lof_map,args.out_path)

    pretty_print("GENES")
    gene_matrix,sample_header = merge_gene_chunks(ordered_genes,lof_dict,args.vcf,args.out_path,args.test,args.force)

    pretty_print("VCF")
    build_vcf(gene_matrix,sample_header,gene_chrom_map,args.out_path)
    
if __name__ == '__main__':

    """
    Builds custom vcf where genotype is probability of carrying at lest one LOF variant.
    """
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
    parser.add_argument('--force',action = 'store_true',help = 'Runs regardless of cache.')

    args=parser.parse_args()
    # logging level
    level = log_levels[args.log]
    logging.basicConfig(level=level,format="%(levelname)s: %(message)s")

    make_sure_path_exists(args.out_path)
    main(args)
