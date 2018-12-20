import argparse
import pandas as pd
from functools import partial
from file_utils import *


###############################
#------MERGE GENES MATRIX-----#
###############################
 
        
@timing_function
def do_chunks(args):
    '''
    Creates the final LOF matrix.
    '''

    print('merging variants into genes...')
    with open(args.g2v,'rb') as i: g2v = pickle.load(i)
    genes = len(list(g2v.keys()))
    print(genes, ' genes present')
    
    # define func that fixes args and accepts only the index of the proc as an argument
    merge_multiproc = partial(merge_gene,args = args)

    pool = multiprocessing.Pool(args.cpus)
    pool.map(merge_multiproc,range(args.cpus))
    pool.close()
    print('done.')

def merge_chunks(args):
    #merge matrix chunks
    print('merging matrix chunks..')
    out_file = os.path.join(args.out_path, args.lof + '_matrix.txt')
    #load samples
    samples = pd.read_csv(args.matrix,sep = args.sep,usecols = ['IID']).values.flatten()    
    with open(out_file,'wt') as o: o.write('GENE\t' + '\t'.join(samples) + '\n')
    print('samples written...')
    cmd = "cat " + os.path.join(args.tmp_path, args.lof +'_' + "matrix_chunk*") + "  >> "+ out_file
    tmp_bash(cmd)
    print('done..')
    
    
def merge_gene(chunk,args):
    '''
    For each chunk of genes, it returns the columns of the variants in the gene (through g2v) and then merges them,keeping 1s were present.
    '''
    matrix_chunk = os.path.join(args.tmp_path, args.lof +'_' + 'matrix_chunk_'+str(chunk) + '.txt')
    if not os.path.isfile(matrix_chunk) or args.force:
        with open(args.g2v,'rb') as i: g2v = pickle.load(i)
        gene_chunk = get_gene_chunk(chunk,args,g2v)
        with open(matrix_chunk,'wt') as o:
            for i,gene in enumerate(gene_chunk):
                header = return_header_variants(args.matrix)
                indexcol =[header.index(element) for element in g2v[gene]]
                #read column --> replace NA with 0s --> convert to INT --> return max val --> bool --> string
                data = pd.read_csv(args.matrix,sep =args.sep,usecols = indexcol).fillna(0).apply(max,axis =1).astype(bool).astype(int).astype(str)
                o.write(gene +'\t' + '\t'.join(data.values) + '\n')
                get_progress(args)
    else:
        print('matrix chunk already generated.')
        
    
def get_gene_chunk(idx,args,g2v):
    """
    Splits the genelist into chunks and returns the chunk matching the index
    """
    
    
    gene_list = np.array(list(g2v.keys()))
    chunk_list = split_array_chunk(gene_list,args.cpus)
    chunk = chunk_list[idx]
    chunk = chunk[:args.test] if args.test else chunk
    return chunk
                        

#############################
#------PROCESS VARIANTS-----#
#############################
def return_lof_variants(args):
    '''
    Creates the following files:
    *_variants.txt : variantlist used later to filter for variants with plink
    *_variants_gene_dict.p : a pickle object containing a variant to gene dict
    '''
    if args.lof == 'hc_lof':
        lof_filter = ["true"]
    elif args.lof == "most_severe":
        lof_filter = ["frameshift_variant","splice_donor_variant","stop_gained","splice_acceptor_variant"]

    args.lof_variants =os.path.join(args.out_path, args.lof + '_variants.txt')
    args.v2g = os.path.join(args.out_path, args.lof + '_variants_gene_dict.p')
    if os.path.isfile(args.lof_variants):
        print('variants already filtered',args.lof_variants)

    else:
        v2g = {}
        print('saving to ', args.lof_variants)
        with gzip.open(args.annotated_file,'rt') as i:header = next(i).strip().split('\t')
        info_iterator = basic_iterator(args.annotated_file,skiprows=1,columns = [header.index(elem) for elem in ['#variant','gene',args.lof]])# iterate through columns that match gene and LOF
        with open(args.lof_variants,'wt') as v:
            for line in info_iterator:
                variant,gene,lof = line
                if (lof in lof_filter):
                    variant =variant.replace(':','_')
                    v2g[variant] = gene
                    v.write(variant+ '\n')
        with open(args.v2g,'wb') as o:
            pickle.dump(v2g,o)
                            
#########################
#------BUILD MATRIX-----#
#########################
def matrix_plink(args):
    '''
    Returns the final list of snps as well as the LOF raw matrix
    '''
   
            
    #write final snplists
    args.snps =  os.path.join(args.tmp_path, args.lof + '.snplist')
    if not os.path.isfile(args.snps) or args.force:
        if args.exclude:
            # merge list of blacklsited variants
            args.exclude_file = os.path.join(args.tmp_path,'excluded_variants.txt')
            print('merging list of variants to exclude...')
            exclude_list = ' '.join(args.exclude)
            cmd = 'cat ' + exclude_list + '| sort | uniq -u > ' + args.exclude_file
            tmp_bash(cmd)
            exclude = ' --exclude '+ pad(args.exclude_file)
        else:
            exclude = ''

        max_maf = '--max-maf ' + pad(str(args.maxMAF)) if args.maxMAF else ''
        cmd = 'plink  -bfile ' + pad(os.path.splitext(args.bed)[0]) +  ' --write-snplist --extract' +pad(args.lof_variants) + pad(exclude) + pad(args.pargs) +  '--threads' + pad(str(args.cpus)) + max_maf +  '--allow-extra-chr --out ' + pad(os.path.splitext(args.snps)[0])
               
        call(shlex.split(cmd))
        args.force = True
    else:
        print('snplist already generated')

        
    # build reordered plink file in order to match the order of the samples passed
    plink_path = os.path.join(args.out_path,'plink')
    make_sure_path_exists(plink_path)
    plink_file = os.path.join(plink_path,args.lof)
    if not os.path.isfile(plink_file + '.bed') or args.force:
        #build fam file out of the args.samples
        args.fam = os.path.join(args.tmp_path, args.lof + '.fam')
        cmd = """awk 'BEGIN{FS=OFS="\t"} {$1 = $1 OFS $1} 1' """ + args.samples + " > " + args.fam
        tmp_bash(cmd)

        cmd = 'plink --bfile ' + pad(os.path.splitext(args.bed)[0]) + '--allow-extra-chr  --extract ' + pad(args.snps) + ' --keep ' + pad(args.fam) +  '--threads ' + pad(str(args.cpus))+' --make-bed  --indiv-sort f' +pad(args.fam) + ' --out ' + pad(plink_file)
        call(shlex.split(cmd))
        cmd = 'plink --bfile ' + pad(plink_file) + '--freq --out ' + pad(plink_file)
        call(shlex.split(cmd))
        args.force = True
    else:
        print('plink file already generated')
        
    #build variant to sample matrix using above mentioned samples
    args.matrix = os.path.join(args.tmp_path, args.lof + '_matrix.raw')
    if not os.path.isfile(args.matrix) or args.force:
              
        remove = '--remove ' + args.remove if args.remove else ''
        cmd = 'plink -bfile ' + pad(plink_file) +  pad(remove) +  '--recode A   --threads' + pad(str(args.cpus)) + '--allow-extra-chr --out' + pad(os.path.splitext(args.matrix)[0])
        print(cmd)
        call(shlex.split(cmd))
        args.force = True
        
    else:
        print('raw matrix already generated')        

    #identify separator of matrix
    with open(args.matrix,'rt') as i: header = i.readline().strip()
    args.sep = identify_separator(header)
    
    

def save_variant_to_gene_dict(args):
    '''
    Reads the final list of variants and dumps a gene to variant dictionary 
    '''
    
    #get variant to gene mapping from full list of variants
    with open(args.v2g,'rb') as i:
        v2g = pickle.load(i)
        
    g2v = dd(list)
    args.g2v = os.path.join(args.out_path, args.lof + '_gene_variants_dict.p')
    if not os.path.isfile(args.g2v) or args.force:
        # read snplist of filtered plink file and keep gene to variant list dictionary
        with open(args.snps,'rt') as i:
            for line in i:
                variant = line.strip()
                gene = v2g[variant]
                g2v[gene].append(variant)

        with open(args.g2v,'wb') as o:
            pickle.dump(g2v,o)
    
if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Deal with lof variants")
    parser.add_argument("--annotated_file", type= file_exists,
                        help="path to annotated file",required = True)
    parser.add_argument("--lof", type= str,
                        help="type of lof filter",required = True,choices = ['hc_lof','most_severe'] )
    parser.add_argument('-o',"--out_path",type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--bed", type= str,help="Path to bed file ",required = True)
    parser.add_argument("--samples", type=file_exists, help =  "Path to list of samples to keep",required = True)

    # not required
    parser.add_argument('--pargs',type = str,default = '',help = 'String with kwargs to pass to all plink calls, based on the desired output')
    parser.add_argument("--maxMAF", type = float, help=  "List of paths to file with list of variants to exclude")
    parser.add_argument("--exclude", type = file_exists, nargs = '*',help=  "List of paths to file with list of variants to exclude")
    parser.add_argument("--remove", type=file_exists, help =  "Path to list of samples to remove")
    parser.add_argument("--cpus", type=int, help =  "Number of cpus to use",default = cpus)
    parser.add_argument("--test", type = int,default = 0, help = "Number of genes to run per cpu, if 0 all genes are run.")
    parser.add_argument('--force',action = 'store_true',help = 'Replaces files by force')
    #parser.add_argument("--cov", type=file_exists, help =  "Covariate file")
    args=parser.parse_args()
    
    make_sure_path_exists(args.out_path)
    args.tmp_path = os.path.join(args.out_path,'tmp')
    make_sure_path_exists(args.tmp_path)
    
    pretty_print('returning' + pad(args.lof) +  'variants')
    return_lof_variants(args)
    
    pretty_print('returning matrix')
    matrix_plink(args)
    save_variant_to_gene_dict(args)
    
    pretty_print('merge genes')
    do_chunks(args)
    merge_chunks(args)
    
