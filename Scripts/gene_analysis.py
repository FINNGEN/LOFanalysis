import numpy as np
import os
import gzip
from collections import defaultdict as dd
import statsmodels.api as sm
from firth_regression import fit_firth
import sys
from scipy.stats import fisher_exact
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
dataPath = rootPath + 'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = rootPath + 'tmp_scripts/'

# REQUIRED FILES
phenoList = np.loadtxt(dataPath + 'pheno-list.txt',usecols = [0],dtype = str,skiprows = 1)
phenoFile = dataPath + 'FINNGEN_PHENOTYPES_DF1_V4_2018_03_27.txt.gz'
eigenvecPath = dataPath + '10pc.eigenvec'


def logistic_regression(iPath,lofString = 'hc_lof',pcData = None,phenoData = None,lofData= None,f = phenoFile,infoFilter = 0):
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

    # here i apply the info score filter
    lofData[lofData > infoFilter] = 1
    y = phenoData
    X = np.c_[lofData,pcData]     
    #logit regression
    logit_model=sm.Logit(y,X)
    result=logit_model.fit()
    #info
    phenoMask = (phenoData >0)
    cases = phenoMask.sum()
    lofCases = lofData[mask].sum()
    
    
    return result,cases,lofCases

def get_lof_data(iPath,gene,lofString = 'hc_lof'):
    with open(iPath + lofString + '_gene_to_filtered_samples.tsv','rt') as i:
        for line in i:
             line = i.readline().strip().split('\t')
             if line[0] == gene:
                 lofData = np.array(line[1:],dtype = float)
                 break
            
    return lofData


def get_pheno_data(iPath,pheno,f = phenoFile,lofString = 'hc_lof'):
    
    data = return_column(pheno = pheno,f = f,dtype = float)
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



    subparsers = parser.add_subparsers(help='help for subcommand',dest ="command")

    # create the parser for the generate_variants command
    parser_fix_samples = subparsers.add_parser('fix-samples', help='fix files in order to match shared samples')
    parser_fix_samples.add_argument("--lof", type= str,
                        help="type of lof filter",required = True )
    parser_fix_samples.add_argument("--oPath", type= str,help="Path to folder where to output",default = ".")
    # create the parser for the "command_2" command
    
    
    args = parser.parse_args()

    if args.command == "fix-samples":
        oPath = (args.oPath + '/' + args.lof +'/').replace('//','/')
        filter_pcs(oPath,args.lof,args.phenoFile,args.pcPath)
        reorder_lof_matrix(oPath,args.lof)
    
    
