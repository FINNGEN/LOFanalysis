import numpy as np
import os
import gzip
from operator import itemgetter
import shutil
from collections import defaultdict
from itertools import product
import pickle
import shlex
from subprocess import Popen, PIPE,call
from file_utils import make_sure_path_exists,return_header_variants,split_array_chunk
import multiprocessing
cpus = multiprocessing.cpu_count()

plink = shutil.which('plink')

rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
dataPath = rootPath + 'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = rootPath + 'tmp_scripts/'
for path in [dataPath,bashPath]:
    make_sure_path_exists(path)


def dd(tp):
    return defaultdict(tp)
def dd_str():
    return defaultdict(str)
lofName = "filtered_lof"
matrixName = "_variantmatrix.tsv"


###############################
#--MERGE VARIANTS INTO GENES--#
###############################

def write_new_matrix(iPath,lofString = 'hc_lof'):

    oFile = iPath + lofString + "_gene_to_sample_lof.tsv"
    if os.path.isfile(oFile):
        print("gene to sample matrix already generated.")

    else:
        g2v = get_variant_to_gene_dict(iPath)
        matrixPath = iPath + lofString + matrixName
        headerVariants = return_header_variants(matrixPath)
        samples =  np.loadtxt(matrixPath,dtype = str,usecols =[0])

        with open(oFile,'wt') as f:
            f.write("\t".join(samples) + '\n')
            for gene in g2v:
                gData = return_gene_columns(gene,iPath,g2v,headerVariants).astype(str)
                gArray = np.concatenate((np.array([gene]),gData))
                assert gArray.shape == samples.shape
                f.write("\t".join(gArray) + '\n')


def do_chunks(iPath,lofString = 'hc_lof'):

    write_genelists(iPath)
    
    params  = list(product(range(cpus),[iPath],[lofString]))
    


def multi_wrapper_func(args):
    multiprocess_func(*args)
    
def multiprocess_func(chunkInt,iPath,lofString):
    
    g2v = get_variant_to_gene_dict(iPath)
    matrixPath = iPath + lofString + matrixName
    headerVariants = return_header_variants(matrixPath)
    chunkPath = iPath + '/gene_chunks/'
    geneList = np.loadtxt(chunkPath + 'gene_chunk_'+str(chunkInt) + '.txt',dtype = str)
    with open(chunkPath + 'matrix_chunk_' + str(chunkInt) + '.tsv','wt') as f:
        for gene in geneList:
             gData = return_gene_columns(gene,iPath,g2v,headerVariants).astype(str)
             gArray = np.concatenate((np.array([gene]),gData))
             f.write("\t".join(gArray) + '\n')
    
def write_genelists(iPath,chunks = cpus):

    chunkPath = iPath + '/gene_chunks/'
    make_sure_path_exists(chunkPath)
    g2v = get_variant_to_gene_dict(iPath)

    geneList = np.array(list(g2v.keys()))
    chunkList = split_array_chunk(geneList,chunks)
    for i,chunk in enumerate(chunkList):
        np.savetxt(chunkPath + 'gene_chunk_'+str(i) + '.txt',chunk,fmt = '%s')
                        


    


    
def return_gene_columns(gene,iPath,g2v,headerVariants,lofString = 'hc_lof'):
    """
    Loops through the header of the matrix file and returns the columns where variants belong to the gene
    """
    geneVariants = g2v[gene]
    matrixPath = iPath +  lofString + '_info_score_matrix.txt'

    geneColumns = [i+1 for i,elem in enumerate(headerVariants) if elem in geneVariants]
    print(gene,geneColumns)

    #import sample data keeping columns of gene
    vData = np.loadtxt(matrixPath,dtype = float,usecols = geneColumns)
    if len(geneColumns) > 1:
        #sum across variants and check if >1
        vData = np.max(vData,axis = 1)
    
    return vData


def get_variant_to_gene_dict(iPath,lofString = 'hc_lof'):
    '''
    Reads the plink snplist and returns a gene to variant dictionary 
    '''
    #get variant to gene mapping from full list of variants
    bFile = iPath +lofString 

    v2g = dd(str)
    with open(dataPath + lofString + '_variants.txt','rt') as i:
        for line in i:
            variant,gene = line.strip().split('\t')
            v2g[variant] = gene

    # read snplist of filtered plink file and keep gene to variant list dictionary
    g2v = dd(list)
    with open(bFile + '.snplist','rt') as i:
        for line in i:
            variant = line.strip()
            gene = v2g[variant]
            g2v[gene].append(variant)
    return g2v




###############################################
#---CONVERTS LOF MATRIX TO INFO_SCORE MATIX---#
###############################################
def write_info_score_matrix(annotatedPath,oPath,lofString,batchPath = dataPath + 'sample_info.txt'):

    '''
    Goes through each line of the matrix(sample data) and updates the 1s to be the Info score for that sample's batch
    '''
    
    oFile = oPath + lofString + '_info_score_matrix.txt'

    if os.path.isfile(oFile):
        print('info score matrix already generated')

    else:
        print('generating info score matrix..')
        matrixPath = oPath + lofString + matrixName

        #stuff required  
        vDict = variant_is_dict(annotatedPath,oPath,lofString)
        s2b = sample_to_batch_ditct(batchPath)
        headerVariants = return_header_variants(matrixPath)
    
        print('looping samples...')
        with open(matrixPath,'rt') as i,open(oFile,'wt') as o:
            next(i) #skip header
            for line in i:
                sample,data = process_line(line,s2b,headerVariants,vDict)
                o.write(sample + ' ' + ' '.join([str(elem) for elem in data]) + '\n')
    

def process_line(line,s2b,headerVariants,vDict):
    '''
    Given a line of the lof_matrix, it returns 0 if not lof and INFO_SCORE of the batch for that variant otherwise
    '''
    line = line.strip().split(' ')
    # get sample info
    sample = line[0]
    batch = s2b[sample]
    #read data from sample
    data = np.array(line[1:],dtype = str)
    data = np.isin(data,['1','2']).astype(float) # boolean if lof or not
    #now we have a 1 if there is lof and 0 elsewhere
    dataMask = np.where(data==1)[0] #index of variant with lof
    # now i create a mini array only for the positions with lof and i multiply the original array, masking it, with the info-score array
    infoArray = np.empty(len(dataMask),dtype = float)
    for i,elem in enumerate(infoArray):
        lofVariant = headerVariants[i]
        infoArray[i] = vDict[lofVariant][batch]
    data[dataMask] *= infoArray

    return sample,data
    



def variant_is_dict(annVariants = annotatedVariants,iPath ='/home/pete/results/hc_lof/',lofString = "hc_lof" ):
    
    '''
    Read the annotated_variants and returns a dict[variant][batch] = INFO_SCORE for teh variants that are in the snplist.
    I can use this dictionary to retreieve the info score for the samples
    '''

    picklePath = dataPath + lofString + '_vDict.p'
    try:
        print('pickling variant/batch/info_score dict..')
        vDict = pickle.load(open(picklePath,'rb'))
    except:
        snplist = iPath + lofString + '.snplist'

        print('data missing, generating..')
        variants = np.loadtxt(snplist,dtype = str)   
        vDict = defaultdict(dd_str)
        with gzip.open(annVariants,'rt') as i:
            #read header
            header = i.readline().strip().split('\t')
            infoPos,lofPos,avgPos,genePos = read_header(header)
            #return position of batches 
            batches = header[infoPos[0]:infoPos[-1]+1]
            #return batches
            batches = [batch.split('INFO_')[1].split('_R1')[0] for batch in batches]
            pickle.dump(batches,open(dataPath + 'ourbatches.p','wb'))
        
            startPos = infoPos[0]
            rangebatches = np.arange(len(batches))
            assert len(batches) == len(infoPos)

            #loop variants
            for line in i:
                line = line.strip().split('\t')
                variant = line[0].replace(':','_')
                if variant in variants:
                    for b in rangebatches:
                        batch = batches[b]
                        vDict[variant][batch] = line[startPos + b]

        pickle.dump(vDict,open(picklePath,'wb'))

    return vDict




def sample_to_batch_ditct(filePath = dataPath + 'sample_info.txt'):
    '''
    Given timo's file maps a sample to a batch. requires a conversion on the fly due to slightly different names between his batch names and ours. Need to pass our batches and use difflib
    '''

    picklePath = dataPath +  'sample_to_batch.p'

    try:
        print('returning sample to batch dict')
        s2b = pickle.load(open(picklePath,'rb'))
        
    except:
        print("data doesn't exist..generating..")
        import difflib                         
        ourbatches = pickle.load(open(dataPath + 'ourbatches.p','rb'))

                             
        s2b = dd(str)
        sampleData = np.loadtxt(filePath,dtype = str,delimiter=':',usecols=[0,-1])
        for entry in sampleData:
            timoBatch,sample = entry
            ourBatch = difflib.get_close_matches(timoBatch,ourbatches)[0]
            s2b[sample] = ourBatch

        pickle.dump(s2b,open(picklePath,'wb'))
    return s2b


#######################
#--GENERATING MATRIX--#
#######################
def generate_matrix(iPath,lofString = 'hc_lof'):
    """
    Returns variant x sample matrix with 1s where variant is present
    """
    iFile = iPath +lofString
    oFile = iFile + matrixName

    if os.path.isfile(oFile):
        print('original lof matrix already generated')

    else:
        cmd = 'plink -bfile '+ iFile +' --recode A --out ' + oFile
        call(shlex.split(cmd))
        #remove unncessary columns
        cmd = "cat " +oFile + ".raw |cut -d ' ' -f-1,7- > " + oFile
        shPath = bashPath +  'filter_lof_matrix.sh'

        with open(shPath,'wt') as o:
            o.write(' #!/bin/bash\n')
            o.write(cmd)

        call(['chmod','+x',shPath])
        call(shPath,shell = True)


def plink_filter(filePath,oPath,geno = 0.9,lofString = "hc_lof"):
    """
    Filter full data for only varianst we need
    """
    snpslist = dataPath + lofString + ".snplist"
    if os.path.isfile(oPath + lofString + ".snplist"):
        print('plink files already generated')

    else:

        make_sure_path_exists(oPath)
        cmd = 'plink -bfile ' + filePath + ' --geno ' + str(geno) + ' --extract ' + snpslist + ' --make-bed -out ' + oPath + lofString
        call(shlex.split(cmd))
        cmd = 'plink -bfile ' + oPath + lofString +  ' --write-snplist --out ' + oPath + lofString
        call(shlex.split(cmd))
    
#-------> here i run wdl


def create_info_file(annotatedFile,lofString = 'hc_lof'):
    '''
    Creates a lof_variants.txt with variants that carry lof along with their genes. 
    '''

    if lofString == 'hc_lof':
        lofFilterList = "true"
    elif lofString == "most_severe":
        lofFilterList = ["frameshift_variant","splice_donor_variant","stop_gained","splice_acceptor_variant"]

    else:
        raise ValueError("invalid lof filter")

    lofPath =dataPath + lofString + '_variants.txt'
    snpsPath =dataPath + lofString + '.snplist'
    if os.path.isfile(lofPath):
        print('variants already filtered')
        return 

    else:
        print('filtering ' + lofString +  ' variants...')
        with gzip.open(annotatedFile,'rt') as i,open(lofPath,'wt') as o,open(snpsPath,'wt') as oo:
            infoPos,lofPos,avgPos,genePos = read_header(i.readline().strip().split('\t'),lofString )

            for line in i:
                line = line.strip().split('\t')
                variant = line[0]
                lof = line[lofPos]
                gene = line[genePos]
                if (lof in lofFilterList):
                    o.write(variant.replace(':','_') + '\t' + gene + '\n')
                    oo.write(variant.replace(':','_')+ '\n')


def read_header(header = None,lofString = "hc_lof"):
    '''
    Reads the first line of the variants.gz file and returns the position of the Info score
    '''
    
    if header == None:
        with gzip.open(annotatedVariants,'rt') as i:
    
            header = i.readline()
            header = header.strip().split('\t')

            
    infoPos = [i for i,elem in enumerate(header) if 'INFO_' in elem]
    lofPos = [i for i,elem in enumerate(header) if  elem == lofString][0]
    avgPos = [i for i,elem in enumerate(header) if 'INFO' == elem][0]
    genePos =  [i for i,elem in enumerate(header) if  elem == "gene"][0]
    return infoPos,lofPos,avgPos,genePos



if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Deal with lof variants")
    parser.add_argument("--annotatedFile", type= str,
                        help="path to annotatedFile",required = False,default =annotatedVariants )

    subparsers = parser.add_subparsers(help='help for subcommand',dest ="command")

    # create the parser for the generate_variants command
    parser_filter = subparsers.add_parser('filter', help='filter the variants')
    parser_filter.add_argument("--lof", type= str,
                        help="type of lof filter",required = True )

    # create the parser for the "command_2" command
    parser_matrix = subparsers.add_parser('generate-matrix', help='help for command_2')
    parser_matrix.add_argument("--plinkPath", type= str,help="Path to plink data (with name of bFile)",required = True)
    parser_matrix.add_argument("--oPath", type= str,help="Path to folder where to output",default = ".")
    parser_matrix.add_argument("--lof", type= str,help="type of lof filter",required = True )
    parser_matrix.add_argument("--geno", type= float,help="genotype call rate for plink",default = 0.9 )


    
    args = parser.parse_args()

    if args.command == "filter":
        create_info_file(args.annotatedFile,args.lof)

    
    if args.command == "generate-matrix":
        oPath = (args.oPath + '/' + args.lof +'/').replace('//','/')
        plink_filter(args.plinkPath,oPath,args.geno,args.lof)
        generate_matrix(oPath,args.lof)
        write_info_score_matrix(args.annotatedFile,oPath,args.lof)
        write_new_matrix(oPath,args.lof)
