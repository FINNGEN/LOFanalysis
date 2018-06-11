import numpy as np
import os
import gzip
from operator import itemgetter
import shutil
from itertools import product
import pickle
import shlex
import sys
from subprocess import Popen, PIPE,call
from file_utils import rootPath,dataPath,annotatedVariants,bashPath
from file_utils import make_sure_path_exists,return_header_variants,split_array_chunk,read_header,get_variant_to_gene_dict,dd_str,sample_to_batch_dict,variant_is_dict
import pandas as pd
import multiprocessing
cpus = multiprocessing.cpu_count()

plink = shutil.which('plink')

lofName = "filtered_lof"
matrixName = "_variantmatrix.tsv"



def matrix_to_bool(iPath,lofString,infoScore=0.9):
    oFile = iPath +'/matrix/' + lofString + '_gene_to_sample_' + str(infoScore) + '.tsv'
    make_sure_path_exists(os.path.dirname(oFile))
    
    mFile = iPath + lofString + "_gene_to_sample.tsv"
    if os.path.isfile(oFile):
        print('matrix already converted')
    else:
        genes =  pd.read_csv(mFile,dtype = str,header = None,sep = '\t',usecols = [0]).values.flatten()
        print('genes loaded.')
        with open(oFile,'wt') as o:
            for i,gene in enumerate(genes):
                sys.stdout.write('processing gene %i %s \r'%(i,gene)),
                sys.stdout.flush()
                line = pd.read_csv(mFile,nrows = 1,skiprows = i,header = None,dtype = str,sep = '\t').values.flatten()[1:].astype(float)
                line = np.where(line < infoScore,0,1).astype(str)
                line = np.concatenate((np.array([gene]),line))
                o.write("\t".join(line) + '\n')

    

###############################
#--MERGE VARIANTS INTO GENES--#
###############################

def write_gene_matrix(iPath,lofString):
    '''
    Merges the chunks
    '''
    oFile = iPath + lofString + "_gene_to_sample.tsv"
    sFile = iPath + lofString + "_samplelist.txt"
    if os.path.isfile(oFile):
        print("gene to sample matrix already generated.")

    else:
        g2v = get_variant_to_gene_dict(iPath,lofString)
        matrixPath = iPath + 'plink_files/'+ lofString + matrixName
        headerVariants = return_header_variants(matrixPath)
        samples =  np.loadtxt(matrixPath,dtype = str,usecols =[0])[1:]
        np.savetxt(sFile,samples,fmt ='%s')
        with open(oFile,'wt') as f:
            chunkPath = iPath + '/gene_chunks/'
            for i in range(cpus):
                sys.stdout.write('\r merging chunk %i'%(i))
                sys.stdout.flush()
                with open(chunkPath + 'matrix_chunk_'+str(i) + '.tsv','rt') as i:
                    for line in i:
                        f.write(line)

def do_chunks(iPath,lofString):
    '''
    Merges the variant file in chunks and outputs a gene_to_sample matrix for each process.
    '''
    write_genelists(iPath,lofString)

    params  = list(product(range(cpus),[iPath],[lofString]))
    pool = multiprocessing.Pool(cpus)
    pool.map(multi_wrapper_func,params)
    pool.close()

def multi_wrapper_func(args):
    multiprocess_func(*args)
        
def multiprocess_func(chunkInt,iPath,lofString):
    '''
    Function that takes a chunk of genes and writes the output matrix
    '''
    chunkPath = iPath + '/gene_chunks/'
    chunkFile = chunkPath + 'matrix_chunk_' + str(chunkInt) + '.tsv'
    if os.path.isfile(chunkFile):
        print('chunk ' + str(chunkInt) + ' already generated')
    else:
        g2v = get_variant_to_gene_dict(iPath,lofString)
        matrixPath = iPath + 'plink_files/'+lofString + matrixName
        headerVariants = return_header_variants(matrixPath)

        geneList = np.loadtxt(chunkPath + 'gene_chunk_'+str(chunkInt) + '.txt',dtype = str)
        with open(chunkFile,'wt') as f:
            for gene in geneList:
                gData = return_gene_columns(gene,iPath,g2v,headerVariants,lofString).astype(str)
                gArray = np.concatenate((np.array([gene]),gData))
                f.write("\t".join(gArray) + '\n')
    
def write_genelists(iPath,lofString,chunks = cpus):
    """
    Splits the genelist into chunks
    """
    chunkPath = iPath + '/gene_chunks/'
    make_sure_path_exists(chunkPath)
    g2v = get_variant_to_gene_dict(iPath,lofString)

    geneList = np.array(list(g2v.keys()))
    chunkList = split_array_chunk(geneList,chunks)
    for i,chunk in enumerate(chunkList):
        np.savetxt(chunkPath + 'gene_chunk_'+str(i) + '.txt',chunk,fmt = '%s')
                        



def return_gene_columns(gene,iPath,g2v,headerVariants,lofString):
    """
    Given a gene it loops through the header of the matrix file and returns the columns where variants belong to the gene. Now i'm adding the feauture to have instead of 1s the highest info score for the batch the sample belongs to.

    # INPUTS
    - gene
    - iPath : inputPath
    - g2v: gene to variant dictionary
    - headerVariants: list of variants to be included.

    # OUTPUTS
    - column(s) of the matrix for the gene
    """

    # path of the matrix file
    matrixPath = iPath + '/plink_files/'+ lofString + matrixName
    # return index of colums if variant belongs to gene
    data = [(i+1,variant) for i,variant in enumerate(headerVariants) if variant in g2v[gene]]
    geneColumns = [elem[0] for elem in data]
    variants = [elem[1] for elem in data]
    print(gene,geneColumns,variants)

    #import sample data keeping columns of gene
    vData = np.loadtxt(matrixPath,dtype = str,usecols = geneColumns,skiprows =1 )
    vData = np.isin(vData,['1','2']) # boolean if lof or not

    # import info score data from same columns
    infoPath = iPath + '/misc/' + lofString + '_info_variant_matrix.tsv'
    iData = np.loadtxt(infoPath,dtype = float,usecols = geneColumns)
    #return maximum info_score value
    res = vData*iData
    if len(geneColumns) > 1:
        res=  np.max(res,axis = 1)
        #sum across variants and check if >1
        #vData = np.sum(vData,axis = 1)

    return res



def write_info_score_matrix_sample(iPath,samplePath,lofString):
    '''
    I build an analogue matrix so that instead of 1s and 0s we have the info_score of the sample

    # INPUTS
    - path to data
    - path to sample to batch data
    - lof

    # OUTPUTS
    - sample to variant matrix with info score of sample's batch as a value
    '''

    print('generating sample to variant info score matrix...')
    oFile = iPath + '/misc/'
    make_sure_path_exists(oFile)
    oFile += lofString + '_info_variant_matrix.tsv'

    if os.path.isfile(oFile):
        print('matrix already generated')

    else:
        matrixPath = iPath + '/plink_files/'+ lofString + matrixName
        samples = np.loadtxt(matrixPath,usecols = [0],dtype = str,skiprows = 1)
    
        s2b = sample_to_batch_dict(samplePath)

        headerVariants = return_header_variants(matrixPath)
    
        vDict =  variant_is_dict(annotatedVariants,iPath,lofString)
        print('data loaded. writing matrix...')
        with open(oFile,'wt') as o:
            for sample in samples:
                batch = s2b[sample]
                line = [sample]
                for variant in headerVariants:
                    line.append(vDict[variant][batch])
                o.write('\t'.join(line) + '\n')                
    return None


#######################
#--GENERATING MATRIX--#
#######################


# THIS PIPELINE GENERATES THE FIRST STEP OF THE PROCESS
# the output is a variant to sample matrix with 1/2/0/NA as entries

def generate_matrix(iPath,lofString = 'hc_lof'):
    """
    Returns variant x sample matrix with 1s where variant is present
    """

    iPath += '/plink_files/'
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
    oPath += '/plink_files/'

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
    It also creates a infoDict, i.e. a dictionary that stores the average info score of the variant
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
            infoDict = dd_str() # dict that store the info score of the variant
            for line in i:
                line = line.strip().split('\t')
                variant = line[0]
                lof = line[lofPos]
                gene = line[genePos]
                avgInfo = line[avgPos]
                if (lof in lofFilterList):
                    variant =variant.replace(':','_') 
                    o.write(variant + '\t' + gene + '\n')
                    oo.write(variant + '\n')
                    infoDict[variant] = float(avgInfo)
                    
        pickle.dump(infoDict,open(dataPath + lofString + '_infoDict.p','wb'))
                    


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Deal with lof variants")
    parser.add_argument("--annotatedFile", type= str,
                        help="path to annotatedFile",required = False,default =annotatedVariants )
    parser.add_argument("--lof", type= str,
                        help="type of lof filter",required = True )

    subparsers = parser.add_subparsers(help='help for subcommand',dest ="command")

    # create the parser for the generate_variants command
    parser_filter = subparsers.add_parser('filter', help='filter the variants')
    # create the parser for the "command_2" command
    parser_matrix = subparsers.add_parser('generate-matrix', help='help for command_2')
    parser_matrix.add_argument("--plinkPath", type= str,help="Path to plink data (with name of bFile)",required = True)
    parser_matrix.add_argument("--oPath", type= str,help="Path to folder where to output",default = ".")
    parser_matrix.add_argument("--geno", type= float,help="genotype call rate for plink",default = 0.9 )
    parser_matrix.add_argument("--samplePath", type = str,help ='path to file with info about sample and batches ',default = dataPath + 'sample_info.txt')

    parser_bool = subparsers.add_parser('bool-matrix', help='help for command_2')
    parser_bool.add_argument("--oPath", type= str,help="Path to folder where to output",default = ".")
    parser_bool.add_argument("--infoScore",type = float,help =' minimum value of info score to convert to bool')

    
    args = parser.parse_args()

    if args.command == "filter":
        create_info_file(args.annotatedFile,args.lof)

    
    if args.command == "generate-matrix":
        oPath = (args.oPath + '/' + args.lof +'/').replace('//','/')
        plink_filter(args.plinkPath,oPath,args.geno,args.lof)
        # build the plink lof matrix
        generate_matrix(oPath,args.lof)
        write_info_score_matrix_sample(oPath,args.samplePath,args.lof)
        #from the plink matrix merge variants into genes
        do_chunks(oPath,args.lof)
        write_gene_matrix(oPath,args.lof)


    if args.command == 'bool-matrix':
        oPath = (args.oPath + '/' + args.lof +'/').replace('//','/')
        matrix_to_bool(args.oPath,args.lof,args.infoScore)
