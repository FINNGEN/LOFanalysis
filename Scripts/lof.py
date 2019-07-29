#!/usr/bin/env python3


from Tools.utils import file_exists,make_sure_path_exists,tmp_bash,pad,return_header,progressBar,mapcount,pretty_print
from file_utils import split_array_chunk,get_progress_test,check_positive
import extract_variants
import argparse,os,multiprocessing,pickle,json
import numpy as np
from itertools import product
import pandas as pd
from collections import defaultdict as dd



def build_gene_matrix(args):

    header_file = os.path.join(args.out_path,args.lof + '_' + str(args.chrom) + '_header.txt')
    #args.samples = np.loadtxt(args.sample_file,dtype = str)
    with open(header_file,'wt') as o:#o.write('GENE\t' + '\t'.join(args.samples) + '\n')
        o.write('GENE')
        with open(args.sample_file,'rt') as i:
            for line in i:
                sample = line.strip()
                o.write('\t' + sample)
        o.write('\n')
    
    args.gene_matrix = os.path.join(args.out_path,args.lof + '_' + str(args.chrom) + '_gene_matrix.tsv')
    if not os.path.isfile(args.gene_matrix) or args.force:
        args.force = True
        pools = multiprocessing.Pool(args.cpus)
        pools.map(gene_multiproc,product([args],range(args.cpus)))
        pools.close()
        matrix_chunk_path = os.path.join(args.out_path,'matrix_chunks/',args.lof + '_CHUNK_gene_matrix.tsv')
        matrix_files = [matrix_chunk_path.replace('CHUNK',str(chunk)) for chunk in range(args.cpus)]
        cmd = f"cat {' '.join(matrix_files)} >> {header_file} && mv {header_file} {args.gene_matrix}"
        tmp_bash(cmd)
        print(' done.')
    else:
        print('gene matrix already calculated')
    assert(mapcount(args.gene_matrix) -1 == len(args.genes))

    
def gene_multiproc(args):
    build_gene_matrix_chunk(*args)
    
def build_gene_matrix_chunk(args,chunk):
    '''
    Merges all variants into their respective genes for a chunk of genes
    '''

    matrix_chunk_path = os.path.join(args.out_path,'matrix_chunks/')
    gene_matrix_file = os.path.join(matrix_chunk_path,args.lof + '_' + str(chunk) + '_gene_matrix.tsv')
    gene_chunk = np.array_split(args.genes,args.cpus)[chunk]
    header = return_header(args.variant_matrix)
    with open(gene_matrix_file,'wt') as out_file:
            for i,gene in enumerate(gene_chunk):
                indexcol =[header.index(element) for element in args.g2v[gene]]
                if args.gp0:
                    data = (1 - pd.read_csv(args.variant_matrix,sep ='\t',usecols = indexcol,dtype = float).prod(axis = 1)).values
                elif args.hard_call:
                    best_variant = 1 - pd.read_csv(args.variant_matrix,sep ='\t',usecols = indexcol,dtype = float).min(axis = 1)
                    data = np.where(best_variant >= args.hard_call,1,0)
                    
                out_file.write(gene +'\t' + '\t'.join(map(str,data)) + '\n')
                get_progress_test(args)


################################
#---VARIANT TO SAMPLE MATRIX---#
################################


def build_raw_matrix(args):
    """
    Function that build the variant to sample matrix for lof variants.
    """

    # if sample file is not provided extract from vcf
    if not args.sample_file:
        args.sample_file = os.path.join(args.out_path,'samples.txt')
        if not os.path.isfile(args.sample_file):
            test_filter = ' | head -n20 ' if args.test else ''
            samples_cmd = f"bcftools query -l {args.vcf} test_filter  > {args.sample_file}"
            tmp_bash(samples_cmd)
    
    args.variant_matrix = os.path.join(args.out_path,args.lof + '_' + str(args.chrom) + '_variant_matrix.tsv')
    if not os.path.isfile(args.variant_matrix) or args.force:
        args.force = True
        build_matrix_chunks(args)
    else:
        print('variant matrix already calculated')

    args.g2v_file = os.path.join(args.variants_path, args.lof + '_gene_variants_dict.p')
    if not os.path.isfile(args.g2v_file) or args.force:
        args.force = True
        save_variant_to_gene_dict(args)
        
    with open(args.g2v_file,'rb') as i: args.g2v = pickle.load(i)
    args.genes = np.array(list(args.g2v.keys()))

    print(len(return_header(args.variant_matrix)), ' lof variants')
    print(len(args.genes), ' lof genes')
    
def build_matrix_chunks(args):
    '''
    Runs a multiproc query for each chunk of lof variants.
    '''
    #path to save chunks
    matrix_chunk_path = os.path.join(args.out_path,'matrix_chunks/')
    make_sure_path_exists(matrix_chunk_path)

    variant_chunk = os.path.join(args.chunk_path,'variant_chunk_CHUNK.txt')

    # sample_filters
    sample_string = f"-S {args.sample_file}" 
    # variant_filters
    variant_filter = f" -i 'ID=@{variant_chunk}"   
    if args.info_score:
        variant_filter += f" & INFO>={args.info_score} "
    
    # AF < 0.5
    pretty_print('AF < 0.5')
    # make sure the variant is rare
    AF_variant_filter =  variant_filter + f" & AF< 0.5' "
        
    matrix_file = os.path.join(matrix_chunk_path,args.lof + '_CHUNK_matrix.tsv')
    chunk_file = os.path.join(args.chunk_path,'position_chunk_CHUNK.tsv')
    #use bcftools to keep positions and info score if rqeuired.
    matrix_cmd = f"bcftools query -R {chunk_file} {sample_string} -f" + " '%ID[\\t%GP{0}]\\n'"  +f"{AF_variant_filter} {args.vcf} > {matrix_file}"   

    pools = multiprocessing.Pool(args.cpus)
    pools.map(multiproc_cmd,product([matrix_cmd],range(args.cpus)))
    pools.close()

    # AF > 0.5
    pretty_print('AF > 0.5')
    AF_variant_filter =  variant_filter + f" & AF > 0.5' "    
    #use bcftools to keep positions and info score if rqeuired.
    matrix_cmd = f"bcftools query -R {chunk_file} {sample_string} -f" + " '%ID[\\t%GP{2}]\\n'"  +f"{AF_variant_filter} {args.vcf} >> {matrix_file}"   

    pools = multiprocessing.Pool(args.cpus)
    pools.map(multiproc_cmd,product([matrix_cmd],range(args.cpus)))
    pools.close()

    pretty_print("SANITY CHECK")
    matrix_cmd = f"cat {matrix_file} |wc -l && cut -f1 {matrix_file} | sort | uniq | wc -l"
    for chunk in range(args.cpus):
        tmp_bash(matrix_cmd.replace('CHUNK',str(chunk)))
        
    print('merging and transposing chunks')        
    matrix_files = [matrix_file.replace('CHUNK',str(i)) for i in range(args.cpus)]
    cmd = f"cat {' '.join(matrix_files)}  | datamash transpose > {args.variant_matrix}"
    tmp_bash(cmd)

    
def multiproc_cmd(args):
    run_cmd(*args)
        
def run_cmd(cmd,chunk):
    '''
    Given a basic bash cmd, it runs them in parallel replacing CHUNK with the chunk number
    '''
    cmd = cmd.replace('CHUNK',str(chunk))
    print(cmd)
    try:
        tmp_bash(cmd)
    except:
        print(chunk,'failed')
    
def save_variant_to_gene_dict(args):
    '''
    Reads the final list of variants and dumps a gene to variant dictionary 
    '''
    
    #get variant to gene mapping from full list of variants
    with open(args.v2g,'rb') as i:v2g = pickle.load(i)
        
    g2v = dd(list)
    final_variants = return_header(args.variant_matrix)
    for variant in final_variants:
        g2v[v2g[variant]].append(variant)
               
    with open(args.g2v_file,'wb') as o:pickle.dump(g2v,o)
    with open(os.path.join(args.out_path,args.lof+'_'+str(args.chrom)+'_gene_variants_dict.txt'),'w') as out_file:
        out_file.write(json.dumps(g2v))

def main(args):
    make_sure_path_exists(args.out_path)
    pretty_print('CHROMOSOME ' +str(args.chrom))
    extract_variants.return_chrom_variants(args)
    build_raw_matrix(args)
    build_gene_matrix(args)

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Deal with lof variants")
    variant_parser = parser.add_mutually_exclusive_group(required = True)
    variant_parser.add_argument("--annotation", type= file_exists,
                        help="path to annotated file")
    variant_parser.add_argument("--lof_variants", type= file_exists,
                        help="path to list of lof variants")

    call_parser = parser.add_mutually_exclusive_group()
    call_parser.add_argument('--hard-call',type = check_positive, help = 'Hard call filter')
    call_parser.add_argument("--gp0", action='store_true', help =  "Works with GP = 0")
    
    parser.add_argument('-o',"--out_path",type = str, help = "folder in which to save the results", required = True)
    parser.add_argument('-c','--chrom',type = int,help = 'chromosome number', required = True)
    
    parser.add_argument("--vcf", type= file_exists,
                        help="path to vcf file",required = True)
    parser.add_argument("--lof", type= str,
                        help="type of lof filter",required = True,choices = ['hc_lof','most_severe'] )
    parser.add_argument("--info_score", type= float, help="Info score filter")
    parser.add_argument('-s',"--sample_file",type = file_exists, help = "list of samples", required = False)
    parser.add_argument('--cpus',type = int,help = 'number of parallel processes to run', default = 2)
    parser.add_argument('--force',action = 'store_true',help = 'Replaces files by force')
    parser.add_argument('--test',action = 'store_true',help = 'Runs test version') 

    args=parser.parse_args()
    args.out_path = os.path.join(args.out_path,str(args.chrom))
    #print(args)

    main(args)

        
