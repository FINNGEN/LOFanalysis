from Tools.utils import basic_iterator,return_header,file_exists,make_sure_path_exists
from file_utils import split_array_chunk
import os,argparse,pickle
from collections import defaultdict as dd

#############################
#------PROCESS VARIANTS-----#
#############################

def split_variants_chunks(args):
    '''
    Splits the lof variants and positions in chunks for multiprocessing with bcftools
    '''
    #read positions
    variants = []
    with open(args.variants,'rt') as i:
        for line in i:variants.append(line)

    args.chunk_path = os.path.join(args.variants_path,'chunks/')
    make_sure_path_exists(args.chunk_path)
    # split positions in chunks and save them to files
    chunks = split_array_chunk(variants,args.cpus)
    for i,chunk in enumerate(chunks):
        pos_chunk = os.path.join(args.chunk_path,'position_chunk_'+str(i) + '.tsv')
        variant_chunk = os.path.join(args.chunk_path,'variant_chunk_'+str(i) + '.txt')
        with open(variant_chunk,'wt') as v,open(pos_chunk,'wt') as p:
            for variant in chunk:
                v.write(variant)
                chrom,pos,*_ = variant.strip().split('_')
                p.write('\t'.join([chrom,pos])+ '\n')


def return_chrom_variants(args):
    '''
    Returns lof variants that belong to the chromosome being run.
    Creates the following files:
    *_variants.txt : variant list used later to filter for variants
    *_positions.txt : position list used later to filter for variants 
    *_variants_gene_dict.p : a pickle object containing a variant to gene dict
    '''
    args.variants_path = os.path.join(args.out_path,'variants/')
    make_sure_path_exists(args.variants_path)
    
    if not args.lof_variants:
        return_lof_variants(args)

    args.positions =os.path.join(args.variants_path, args.lof + '_' + str(args.chrom) + '_positions.txt')
    args.variants =os.path.join(args.variants_path, args.lof + '_' + str(args.chrom) + '_variants.txt')
    args.v2g = os.path.join(args.variants_path, args.lof + '_'+str(args.chrom) + '_variants_gene_dict.p')
    if not os.path.isfile(args.positions) or not os.path.isfile(args.variants) or not os.path.isfile(args.v2g) or args.force:
        args.force = True
        v2g = {}
        startChrom = 'chr' + str(args.chrom) if args.chrom != 23 else 'chrX'
        iterator = basic_iterator(args.lof_variants,separator = '\t')
        with open(args.positions,'wt') as p,open(args.variants,'wt') as v:                 
            for variant,gene in iterator:
                chrom,pos,ref,alt = variant.split('_')
                if chrom == startChrom:
                    v.write('_'.join((chrom,pos,ref,alt)) +'\n')
                    p.write('\t'.join((chrom,pos)) + '\n')
                    v2g[variant] = gene
        with open(args.v2g,'wb') as o:
            pickle.dump(v2g,o)

    split_variants_chunks(args)
    print("variant chunks created")

    
def return_lof_variants(args):
    '''
    Extracts lof variants from annotation file
    '''
    if args.lof == 'hc_lof':
        lof_filter = ["true"]
    elif args.lof == "most_severe":
        lof_filter = ["frameshift_variant","splice_donor_variant","stop_gained","splice_acceptor_variant"]

    args.lof_variants =os.path.join(args.variants_path, args.lof + '_variants.txt')
    print('saving lof variants to  ', args.lof_variants)
    header = return_header(args.annotation)
    info_iterator = basic_iterator(args.annotation,skiprows=1,columns = [header.index(elem) for elem in ['#variant','gene',args.lof]])

    with open(args.lof_variants,'wt') as o:
        for variant,gene,lof in info_iterator:
            if (lof in lof_filter):
                variant =variant.replace(':','_')
                o.write(variant + '\t' + gene + '\n')

    
if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Deal with lof variants")
    variant_parser = parser.add_mutually_exclusive_group(required = True)
    variant_parser.add_argument("--annotation", type= file_exists,
                        help="path to annotated file")
    variant_parser.add_argument("--lof_variants", type= file_exists,
                        help="path to list of lof variants")
    
    parser.add_argument("--lof", type= str,
                        help="type of lof filter",required = True,choices = ['hc_lof','most_severe'] )
    parser.add_argument('-o',"--out_path",type = str, help = "folder in which to save the results", required = True)
    parser.add_argument('-c','--chrom',type = int,help = 'chromosome number', required = True)
    parser.add_argument('--cpus',type = int,help = 'number of parallel processes to run', default = 2)
    parser.add_argument('--force',action = 'store_true',help = 'Replaces files by force') 

    args=parser.parse_args()
    print(args)
    make_sure_path_exists(args.out_path)
    return_chrom_variants(args)
