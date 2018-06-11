import numpy as np
import os
from collections import defaultdict as dd
import sys
import shlex
from itertools import product
from scipy.stats import fisher_exact
from file_utils import make_sure_path_exists,return_column, get_variant_to_gene_dict,f_test
from file_utils import dataPath,phenoFile
import pandas as pd
import multiprocessing
cpus = multiprocessing.cpu_count()


# REQUIRED FILES
phenoList = np.loadtxt(dataPath + 'pheno-list.txt',usecols = [0],dtype = str,skiprows = 1)

#####################
#--GENE  MULTIPROC--#
#####################


def fisher_gene(iPath,lofString='hc_lof',f = phenoFile,proc = cpus,test = True,infoScore = 0.9):
    '''
    In order to speed up reading from disk, I multiproc the pheno data within a gene.
    '''
    matrix = iPath + lofString + "_gene_to_sample.tsv"
    geneList = np.genfromtxt(matrix,usecols = [0],dtype = str,delimiter = '\t')
    print(len(geneList))
    gList = geneList if test is False else geneList[:10]  
    print(len(gList),' genes')

    pList = phenoList if test is False else phenoList[:proc]
    for gene in gList:
        gene_proc(iPath,pList,lofString,gene,f,proc,infoScore)

def gene_proc(iPath,phenoList,lofString='hc_lof',gene = 'TTLL10',f = phenoFile,proc = cpus,infoScore = 0.9):
    '''
    For a single gene, it multiprocs over the pheno data
    '''

    # where to write_Results
    oPath = iPath + '/gene_fisher/'
    make_sure_path_exists(oPath)
    oFile = oPath + lofString + '_' + gene + '_' + str(infoScore)+  '_results.txt'

    # load lofData
    print(gene)
    lofData = get_lof_data(iPath,gene,lofString,infoScore)
 
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
    return fisher_gene(*args)
def fisher_gene(iPath,lofData,pheno,lofString = 'hc_lof',f = phenoFile):
    '''
    Function that is ultimately passed to the multiprocessing pool. It performs a fisher test given lofData and a phenotype.
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




###############
#--READ DATA--#
###############

def get_lof_data(iPath,gene,lofString = 'hc_lof',infoScore = 0.9):
    '''
    Given a gene, it reads the row of the matrix
    '''
    matrix = iPath +'/matrix/' + lofString + '_gene_to_sample_' + str(infoScore) + '_shared.tsv'

    with open(matrix,'rt',os.O_NONBLOCK) as i:
        for line in i:
             line = line.strip().split('\t')
             if line[0] == gene:
                 lofData = np.array(line[1:],dtype = int)
                 break
            
    return lofData

def get_pheno_data(iPath,pheno,f = phenoFile,lofString = 'hc_lof'):
    '''
    Given a pheno, it returns the pheno data that matches the lofSamples
    '''
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



#####################################
#--FIX FILES TO ORDER SAMPLE DATA---#
#####################################


def reorder_lof_matrix(iPath,lofString = 'hc_lof',infoScore = 0.9):
    '''
    Shuffles the column of the gene_lof_matrix so that the ordering of shared samples is the same as in the eigenvec file
    '''
    
    
    iMatrix = iPath +'/matrix/' + lofString + '_gene_to_sample_' + str(infoScore) + '.tsv'
    oMatrix = iPath +'/matrix/' + lofString + '_gene_to_sample_' + str(infoScore) + '_shared.tsv'

    if os.path.isfile(oMatrix):
        print('matrix already reordered')
    else:
        samples = get_shared_samples(iPath,lofString)
        lofSamples = return_lof_samples(iPath,lofString)
        with open(iMatrix,'rt') as i,open(oMatrix ,'wt') as o:
            for line in i:
                line = line.strip().split('\t')
                # i create a dict that stores the info for each sample so i can then rearrange them in the pc sample order
                newLine = process_line(line,lofSamples,samples)
                o.write(newLine + '\n')


def process_line(line,lofSamples,samples):
    geneDict = dd()
    #gene info
    gene = line[0]
    sys.stdout.write('processing gene %s \r'%(gene)),
    sys.stdout.flush()
    # store sample data in dict
    data = line[1:]
    assert len(data) == len(lofSamples)
    # store sample Data
    for j,sample in enumerate(lofSamples):
        geneDict[sample] = data[j]
    # keep only lines of samples that belong to the shared sample list
    newLine = gene +'\t'+ '\t'.join([geneDict[tmpSample] for tmpSample in samples]) # based on lof data order!     
    return newLine

###############################
#----RETURN SHARED SAMPLES----#
###############################

#Here I save an np.array with the shared samples based on the order in which they appear in the matrix

def get_shared_samples(iPath,lofString = 'hc_lof',f = phenoFile):
    '''
    Returns and saves the samples shared across all files and it returns them in the order of the lofMatrix.
    After the first reordering, the function will always load the same order of samples
    '''
    sharedPath = iPath + '/misc/' + lofString + '_shared_samples.txt'

    if os.path.isfile(sharedPath):
        print('loading data from ' + sharedPath)
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
        finalSamples = np.array([s for s in lofSamples if s in samples])
        np.savetxt(sharedPath,finalSamples,fmt ='%s')
    return finalSamples


def return_lof_samples(iPath,lofString = 'hc_lof'):
    
    sFile = iPath + lofString + "_samplelist.txt"
    samples = pd.read_csv(sFile,header = None).values.flatten()
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
    parser_fix_samples.add_argument('--infoScore',type = float,default = 0.9, help = 'info score filter value')

    # create the parser for the fisher command
    parser_fisher = subparsers.add_parser('fisher', help='do fisher analysis')  
    parser_fisher.add_argument("--cpus",type = int, help = 'Number of cores to use', default = cpus)
    parser_fisher.add_argument('--test',action = 'store_true',help = 'Flag to run small chunks')
    parser_fisher.add_argument('--infoScore',type = float,default = 0.9, help = 'info score filter value')

    args = parser.parse_args()
    oPath = (args.oPath + '/' + args.lof +'/').replace('//','/')

    if args.command == "fix-samples":
        reorder_lof_matrix(oPath,args.lof,args.infoScore)
    
    
    if args.command == "fisher":
        fisher_gene(oPath,args.lof,args.phenoFile,args.cpus,args.test,args.infoScore)
