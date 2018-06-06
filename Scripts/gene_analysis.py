import numpy as np
import os
from collections import defaultdict as dd
import sys
import shlex
from itertools import product
from scipy.stats import fisher_exact
from file_utils import make_sure_path_exists,return_column
from file_utils import dataPath
from filter_variants import get_variant_to_gene_dict
import pandas as pd
import multiprocessing
cpus = multiprocessing.cpu_count()



# REQUIRED FILES
phenoList = np.loadtxt(dataPath + 'pheno-list.txt',usecols = [0],dtype = str,skiprows = 1)
phenoFile = dataPath + 'FINNGEN_PHENOTYPES_DF1_V4_2018_03_27.txt.gz'

#####################
#--GENE  MULTIPROC--#
#####################


def logit_gene(iPath,lofString='hc_lof',f = phenoFile,proc = cpus,test = True,infoFilter = 0.9):
    '''
    In order to seepd up reading from disk, I multiproc the pheno data within a gene.
    '''
 
    geneList = get_info_score_gene_list(iPath,lofString,infoFilter)
    print(len(geneList))
    gList = geneList if test is False else geneList[:10]  
    print(len(gList),' genes')

    pList = phenoList if test is False else phenoList[:proc]
    for gene in gList:
        gene_proc(iPath,pList,lofString,gene,f,proc)

def gene_proc(iPath,phenoList,lofString='hc_lof',gene = 'TTLL10',f = phenoFile,proc = cpus):
    '''
    For a single gene, multiproc the phenodata
    '''

    # where to write_Results
    oPath = iPath + '/gene_fits/'
    make_sure_path_exists(oPath)
    oFile = oPath + lofString + '_' + gene +  '_gene_results.txt'

    # load lofData
    print(gene)
    lofData = get_lof_data(iPath,gene,lofString)
 
    # list of params to loop
    params  = list(product([iPath],[lofData],phenoList,[lofString],[f]))
    pool = multiprocessing.Pool(proc)
    results = pool.map(gene_wrapper,params)
    pool.close()
    
    with open(oFile,'wt') as o:
        o.write('\t'.join(shlex.split('gene lof_cases lof_controls no_lof_cases no_lof_controls odds ratio fischer_pval ')) + '\n')

        for entry in results:
            pheno,f_results,table = entry
            #print(pheno,i,gene)
            o.write(pheno + '\t')
            # write counts of lof/no_lof
            countString =  '\t'.join([str(elem) for elem in table.flatten()])
            o.write(countString + '\t')
            o.write( '\t'.join([str(elem) for elem in f_results]))
            o.write('\n')
            
        
    return None


def gene_wrapper(args):
    return logistic_gene(*args)
def logistic_gene(iPath,lofData,pheno,lofString = 'hc_lof',f = phenoFile):
    '''
    Function that is ultimately passed to the multiprocessing pool. It loops through all genes given a phenotype. With test  it only works with a small chunk of genes
    '''
 
    phenoDataPath = iPath + '/pheno_data/'
    make_sure_path_exists(phenoDataPath)
    phenoSave = phenoDataPath + lofString + '_' + pheno  + '_phenodata.txt'
    try:
        phenoData = pd.read_csv(phenoSave,header = None).values.flatten()

    except:
        phenoData = get_pheno_data(iPath,pheno,f,lofString)
        np.savetxt(phenoSave,phenoData,fmt = '%i')
        
    f_results,table = f_test(phenoData,lofData)
    return pheno,f_results,table
          


def f_test(phenoData,lofData):
    # get lof counts for cases
    phenoMask = (phenoData >0)
    cases = phenoMask.sum()
    lofCases = int(lofData[phenoMask].sum())
    nolofCases = cases - lofCases
    # get lof counts for controls
    phenoMask = (phenoData == 0)
    controls = phenoMask.sum()
    lofControls = int(lofData[phenoMask].sum())
    nolofControls = controls - lofControls
    # do fischer test
    table = np.empty((2,2),dtype = int)
    table[0] = [lofCases,lofControls]
    table[1] = [nolofCases,nolofControls]
    f_results = fisher_exact(table)

    return f_results,table

def get_lof_data(iPath,gene,lofString = 'hc_lof'):
    with open(iPath + lofString + '_gene_to_filtered_samples.tsv','rt',os.O_NONBLOCK) as i:
        for line in i:
             line = line.strip().split('\t')
             if line[0] == gene:
                 lofData = np.array(line[1:],dtype = int)
                 break
            
    return lofData

def get_pheno_data(iPath,pheno,f = phenoFile,lofString = 'hc_lof'):

    # read data and fix nans by setting all non 0 values
    data = return_column(pheno = pheno,f = f,dtype = float)
    mask = (data >0)
    data[mask] = 1
    data[~mask] = 0
    
    phenoSamples= return_column(f =f,dtype =str)
    phenoDict = dd()
    for i,entry in enumerate(data):
        phenoDict[phenoSamples[i]] = entry
        
    samples = get_shared_samples(iPath,lofString)
    phenoData = np.empty_like(samples,dtype = int)
    for i,sample in enumerate(samples):
        phenoData[i] = int(phenoDict[sample])
    return phenoData

def get_info_score_gene_list(iPath,lofString = 'hc_lof',infoFilter = 0.9):
    '''
    Given a lof and infoFilter it returns the list of genes for that lof that all variants above infoFilter
    '''
    import pickle
    g2v = get_variant_to_gene_dict(iPath,lofString)
    infoDict = pickle.load(open(dataPath + lofString + '_infoDict.p','rb'))
    geneList = [gene for gene in g2v if np.min([infoDict[variant] for variant in g2v[gene]]) > infoFilter]
    return geneList


#####################################
#--FIX FILES TO ORDER SAMPLE DATA---#
#####################################


def reorder_lof_matrix(iPath,lofString = 'hc_lof'):
    '''
    Shuffles the column of the gene_lof_matrix so that the ordering of shared samples is the same as in the eigenvec file
    '''
    
    iMatrix = iPath + lofString + "_gene_to_sample.tsv"
    oMatrix = iPath + lofString + "_gene_to_filtered_samples.tsv"

    if os.path.isfile(oMatrix):
        print('matrix already reordered')
    else:
        samples = get_shared_samples(iPath,lofString)
        lofSamples = return_lof_samples(iPath,lofString)
        with open(iMatrix,'rt') as i,open(oMatrix ,'wt') as o:
            for line in i:
                line = line.strip().split('\t')
                # i create a dict that stores the info for each sample so i can then rearrange them in the pc sample order
                newline = process_line(line,lofSamples,samples)
                o.write(newLine + '\n')


def process_line(line,lofSsamples,samples):
    geneDict = dd()
    #gene info
    gene = line[0]
    sys.stdout.write('processing gene %s \r'%(gene)),
    sys.stdout.flush()
    # sample data
    data = line[1:]
    assert len(data) == len(lofSamples)
    # store sample Data
    for j,sample in enumerate(lofSamples):
        geneDict[sample] = data[j]
        # write new line
    newLine = gene +'\t'
    newline += '\t'.join([geneDict[tmpSample] for tmpSample in samples]) # based on pc data order!
    return newline

###############################
#----RETURN SHARED SAMPLES----#
###############################

#Here I save an np.array with the shared samples based on the order in which they appear in the matrix

def get_shared_samples(iPath,lofString = 'hc_lof',f = phenoFile):
    '''
    Returns and saves the samples shared across all files and it returns them in the order of the lofMatrix.
    After the first reordering, the function will always load the same order of samples
    '''
    sharedPath = iPath +lofString + '_shared_samples.txt'
    if os.path.isfile(sharedPath):
        finalSamples = pd.read_csv(sharedPath,header = None,dtype = str).values.flatten()

    else:
        print('importing all samples..')
        lofSamples = set(return_lof_samples(iPath,lofString))
        phenoSamples = set(return_column(f =f,dtype = str))
        print('done')
        # merging of samples to create unique set
        samples = list(lofSamples.intersection(phenoSamples))
        print('set of samples calculated')
        # keep lof samples in order
        finalSamples = [s for s in lofSamples if s in samples]
        np.savetxt(sharedPath,finalSamples,fmt ='%s')
    return finalSamples


def return_lof_samples(iPath,lofString = 'hc_lof'):
    matrixPath =  iPath + lofString + "_gene_to_sample.tsv"
    with open(matrixPath,'rt') as i:
        samples =i.readline().strip().split('\t')[1:]
    return np.array(samples)





if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Analyze lof data")
    parser.add_argument("--phenoFile", type= str,
                        help="path to phenotype file",required = False,default = phenoFile)
    parser.add_argument("--lof", type= str,
                        help="type of lof filter",required = True )
    parser.add_argument("--oPath", type= str,help="Path to folder where to output",default = ".")

    subparsers = parser.add_subparsers(help='help for subcommand',dest ="command")

    # create the parser for fixing samples
    parser_fix_samples = subparsers.add_parser('fix-samples', help='fix files in order to match shared samples')

    # create the parser for the logit command
    parser_logit = subparsers.add_parser('logit', help='do logit analysis')  
    parser_logit.add_argument("--cpus",type = int, help = 'Number of cores to use', default = cpus)
    parser_logit.add_argument('--test',action = 'store_true',help = 'Flag to run small chunks')
    parser_logit.add_argument('--infoFilter',type = float,default = 0.9, help = 'info score filter value')
    args = parser.parse_args()
    oPath = (args.oPath + '/' + args.lof +'/').replace('//','/')

    if args.command == "fix-samples":
        reorder_lof_matrix(oPath,args.lof)
    
    
    if args.command == "logit":
#        multiproc_logit(oPath,args.lof,args.phenoFile,args.cpus,args.test,args.infoFilter)

        logit_gene(oPath,args.lof,args.phenoFile,args.cpus,args.test,args.infoFilter)
