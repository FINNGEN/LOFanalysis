import numpy as np
import os
import gzip
from collections import defaultdict as dd
import statsmodels.api as sm
from firth_regression import fit_firth
import sys
import shlex
from itertools import product
from scipy.stats import fisher_exact
from file_utils import make_sure_path_exists
from filter_variants import get_variant_to_gene_dict
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
import multiprocessing
cpus = multiprocessing.cpu_count()

rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
dataPath = rootPath + 'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = rootPath + 'tmp_scripts/'

# REQUIRED FILES
phenoList = np.loadtxt(dataPath + 'pheno-list.txt',usecols = [0],dtype = str,skiprows = 1)
phenoFile = dataPath + 'FINNGEN_PHENOTYPES_DF1_V4_2018_03_27.txt.gz'
eigenvecPath = dataPath + '10pc.eigenvec'



def multiproc_logit(iPath,lofString='hc_lof',f = phenoFile,proc = cpus,test = True):

    pList = phenoList if test is False else phenoList[:proc]
    pList = phenoList
    print(len(pList))
    params  = list(product([iPath],pList,[lofString],[f],[test]))
    pool = multiprocessing.Pool(proc)
    pool.map(logit_wrapper,params)
    pool.close()
 
def logit_wrapper(args):
    logistic_pheno(*args)
def logistic_pheno(iPath,pheno,lofString = 'hc_lof',f = phenoFile,test = True):

    oPath = iPath + '/fits/'
    make_sure_path_exists(oPath)
    oFile = oPath + lofString + '_' + pheno + '_pheno_results.txt'
 
    print(pheno)
    phenoData = get_pheno_data(iPath,pheno,f,lofString)
    print('phenoData imported')
    pcPath = iPath + lofString + '_pcs.txt'
    pcData = np.loadtxt(pcPath,dtype = float,usecols = range(1,11))
    print('pcData imported.')
    
    with open(oFile,'wt') as o:
        o.write('\t'.join(shlex.split('gene lof_cases lof_controls no_lof_cases no_lof_controls logit_coeff_gene logit_pval_gene logit_coeff_pc1 logit_pval_pc1 logit_coeff_pc2 logit_pval_pc2 fischer_oddsratio fischer_pval ')) + '\n')
        geneList = get_gene_list(iPath,lofString)
        if test is True:
            geneList = geneList[:100]
        print(len(geneList))
        for i,gene in enumerate(geneList):
            print(pheno,i,gene)
            o.write(gene + '\t')
            lofData = get_lof_data(iPath,gene,lofString)
            logit_results,f_results,table = logistic_regression(iPath,lofString,pcData,phenoData,lofData,f)
            # write counts of lof/no_lof
            countString =  '\t'.join([str(elem) for elem in table.flatten()])
            o.write(countString + '\t')
            
            # write logit_results

            try:
                 params = logit_results.params
                 pvalues = logit_results.pvalues
                 #add columns
                 res = np.column_stack((params,pvalues))
                 # flatten so first two elemts are from lof, next 2 pc1 etc.
                 resArray = res.flatten()[:6]
            except:
                resArray = ['NA']*6
            
            oString =  '\t'.join([str(elem) for elem in resArray])
            o.write( oString + '\t')
            o.write( '\t'.join([str(elem) for elem in f_results]))
            o.write('\n')
                                
            
    
    return None
def logistic_regression(iPath,lofString = 'hc_lof',pcData = None,phenoData = None,lofData= None,f = phenoFile):
    '''
    Returns the logistic regression for the lof + pcs vs pheno.
    phenoDict and lofDictmap samples to their respective value. I need it in order to build arrays that in sync with the pc data
    '''

    if lofData is None:
        print('lofData missing, creating..')
        gene = 'TTLL10'
        print(gene)
        lofData = get_lof_data(iPath,gene,lofString)
        print('done.')
        
    if phenoData is None:
        print('phenoData missing, creating..')
        pheno = phenoList[0]
        print(pheno)
        phenoData = get_pheno_data(iPath,pheno,f,lofString)
        print('done.')


    #now i upload the pc data,along with the samples
    if pcData is None:
        print('importing pcData...')
        pcPath = iPath + lofString + '_pcs.txt'
        pcData = np.loadtxt(pcPath,dtype = float,usecols = range(1,11))

   
    y = phenoData
    X = np.c_[lofData,pcData]     

    #logit regression
    logit_results = None
    # don't run test if there are no lof cases
    if lofData.sum() == 0:
        pass
    else:
        try:
            logit_model=sm.Logit(y,X)
            logit_results=logit_model.fit(full_output = 0,disp = 0)
        except:
            pass
    #fischer test
    samples = len(phenoData)
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

    
    return logit_results,f_results,table

def get_lof_data(iPath,gene,lofString = 'hc_lof'):
    with open(iPath + lofString + '_gene_to_filtered_samples.tsv','rt') as i:
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
    import pickle
    g2v = get_variant_to_gene_dict(iPath,lofString = 'hc_lof')
    infoDict = pickle.load(open(dataPath + lofString + '_infoDict.p','rb'))
    geneList = [gene for gene in g2v if np.max([infoDict[variant] for variant in g2v[gene]]) > infoFilter]
    return geneList
    
def get_gene_list(iPath,lofString = 'hc_lof'):

    oFile = iPath + lofString + "_genelist.txt"
    try:
        geneList = np.loadtxt(oFile,dtype = str)
    except:
        oMatrix = iPath + lofString + "_gene_to_filtered_samples.tsv"
        geneList = np.loadtxt(oMatrix,usecols = [0], skiprows = 1,dtype = str)
        np.savetxt(oFile,geneList,fmt = '%s')
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
                # i create a dict that stores the info for each sample so i can then rearrange them in the pc sample order
                geneDict = dd()
                line = line.strip().split('\t')
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
                newLine += '\t'.join([geneDict[tmpSample] for tmpSample in samples]) # based on pc data order!
                o.write(newLine + '\n')


        
def return_pc_samples(iPath,lofString = 'hc_lof'):

    pcPath = iPath + lofString + '_pcs.txt'
    pcSamples = np.loadtxt(pcPath,dtype = str,usecols = [0])
    return pcSamples
    
    
def filter_pcs(iPath,lofString='hc_lof',f = phenoFile,pcPath = eigenvecPath):
    '''
    Filters the eigenvec file to keep only samples that are shared across all files
    '''

    pcFile = iPath + lofString + '_pcs.txt'
    print('filtering pcs to shared samples --> ' + pcFile)
    if os.path.isfile(pcFile):
        print('principal components already filtered')
    else:
        samples = get_shared_samples(iPath,lofString,f, pcPath)
        print('samples loaded.')
        with open(pcPath,'rt') as i,open(pcFile,'wt') as o:
            for line in i:
                sample = line.strip().split(' ')[0]
                if sample in samples:
                    o.write(line)




###############################
#----RETURN SHARED SAMPLES----#
###############################

#Here I save an np.array with the shared samples based on the order in which they appear in the eigenvector data

def get_shared_samples(iPath,lofString = 'hc_lof',f = phenoFile,pcPath = eigenvecPath):
    '''
    Returns and saves the samples shared across all files and it returns them in the order of the pc file.
    After the first reordering, the function will always load the same order of samples
    '''
    sharedPath = iPath +lofString + '_shared_samples.txt'
    if os.path.isfile(sharedPath):
        finalSamples = np.loadtxt(sharedPath,dtype = str)
    else:
        print('importing all samples..')
       

        def return_pc_original_samples(pcPath = eigenvecPath):
            return np.loadtxt(pcPath,dtype = str,usecols = [0])
       
        pcSamples = return_pc_original_samples(pcPath)
        lofSamples = return_lof_samples(iPath,lofString)
        phenoSamples = return_column(f =f,dtype = str)
        print('done')
        sampleLists = [pcSamples,lofSamples,phenoSamples]
        samples = set(sampleLists[0])
        for sampleList in sampleLists[1:]:
            samples = samples.intersection(sampleList)

        samples = np.array(list(samples))
        print('set of samples calculated')
        finalSamples = [s for s in pcSamples if s in samples]
        np.savetxt(sharedPath,finalSamples,fmt ='%s')
    return finalSamples


def return_lof_samples(iPath,lofString = 'hc_lof'):
    matrixPath =  iPath + lofString + "_gene_to_sample.tsv"
    with open(matrixPath,'rt') as i:
        samples =i.readline().strip().split('\t')[1:]
    return np.array(samples)

##############################
#--PARSE THE PHENOTYPE FILE--#
##############################

def return_column(pheno = 'FINNGENID',f = phenoFile,dtype = 'f8'):

    header = return_header(f =f )
    for i,elem in enumerate(header):
        if str(elem) == pheno:
            phenocol = i
    idcol = 0
    if f.split('.')[-1] == 'txt':
        i = open(f,'rb')
    elif f.split('.')[-1] == 'gz':
        i = gzip.open(f,'rb')

    column = np.genfromtxt(i,usecols = (phenocol,),delimiter = ('\t'),skip_header=1,dtype = dtype)
    i.close()
    return column




def return_header(f = phenoFile):
    '''
    Reads the header of the pheno file
    '''
    if f.split('.')[-1] == 'txt':
        i = open(f,'rt')

    elif f.split('.')[-1] == 'gz':
        i = gzip.open(f,'rt')
    
    header = i.readline()
    header = header.strip().split('\t')
    i.close()
    return header



if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Analyze lof data")
    parser.add_argument("--phenoFile", type= str,
                        help="path to phenotype file",required = False,default = phenoFile)
    parser.add_argument("--pcPath", type= str,
                        help="path to eigenvec file",required = False,default = eigenvecPath)
    parser.add_argument("--lof", type= str,
                        help="type of lof filter",required = True )
    parser.add_argument("--oPath", type= str,help="Path to folder where to output",default = ".")

    subparsers = parser.add_subparsers(help='help for subcommand',dest ="command")

    # create the parser for the generate_variants command
    parser_fix_samples = subparsers.add_parser('fix-samples', help='fix files in order to match shared samples')



    # create the parser for the generate_variants command
    parser_logit = subparsers.add_parser('logit', help='do logit analysis')
   
    parser_logit.add_argument("--cpus",type = int, help = 'Number of cores to use', default = cpus)
    parser_logit.add_argument('--test',action = 'store_true',help = 'Flag to run small chunks')

        # create the parser for the generate_variants command
    parser_pheno = subparsers.add_parser('write-pheno', help='write phenoData')
    parser_pheno.add_argument("--cpus",type = int, help = 'Number of cores to use', default = cpus)
    parser_pheno.add_argument('--test',action = 'store_true',help = 'only runs one pheno per cpu passed')



    args = parser.parse_args()
    oPath = (args.oPath + '/' + args.lof +'/').replace('//','/')

    if args.command == "fix-samples":
        filter_pcs(oPath,args.lof,args.phenoFile,args.pcPath)
        reorder_lof_matrix(oPath,args.lof)
    
    
    if args.command == "logit":

        multiproc_logit(oPath,args.lof,args.phenoFile,args.cpus,args.test)

    if args.command == "write-pheno":

        multiproc_write_pheno(oPath,args.lof,args.phenoFile,args.cpus,args.test)
