import numpy as np
import os
import gzip
from collections import defaultdict as dd
rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
dataPath = rootPath + 'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = rootPath + 'tmp_scripts/'

# REQUIRED FILES
phenoList = np.loadtxt(dataPath + 'pheno-list.txt',usecols = [0],dtype = str,skiprows = 1)
phenoFile = dataPath + 'FINNGEN_PHENOTYPES_DF1_V4_2018_03_27.txt.gz'
phenoFile = dataPath + 'FINNGEN_PHENOTYPES_DF1_2018_03_01.txt'
eigenvecPath = dataPath + '10pc.eigenvec'





def logistic_regression(iPath,lofString = 'hc_lof',phenoDict = None,lofDict= None,f = phenoFile):
    '''
    Returns the logistic regression for the lof + pcs vs pheno.
    phenoDict and lofDictmap samples to their respective value. I need it in order to build arrays that in sync with the pc data
    '''

    if lofDict is None:
        print('lofDict missing, creating...')
        lofDict= dd()
        gene = 'TTLL10'
        print(gene)
        lofSamples = return_lof_samples(iPath,lofString)
        with open(iPath + lofString + '_gene_to_sample.tsv') as i:
            next(i)
            line = i.readline().strip().split('\t')
            assert line[0] == gene
            data = np.array(line[1:],dtype = float)
            assert data.shape == lofSamples.shape
        for i,entry in enumerate(data):
            lofDict[lofSamples[i]] = entry
        print('done.')
        
    if phenoDict is None:
        print('phenoDict missing, creating..')
        pheno = phenoList[0]
        print(pheno)
        data = return_column(pheno = pheno,f = f,dtype = float)
        samples= return_column(f =f,dtype =str)
        assert data.shape == samples.shape
        phenoDict = dd()
        for i,entry in enumerate(data):
            phenoDict[samples[i]] = entry

    #now i upload the pc data,along with the samples
    pcPath = iPath + lofString + '_pcs.txt'
    pcSamples =return_pc_samples(pcPath)
    pcData = np.loadtxt(pcPath,dtype = float,usecols = range(1,11))

    
    lofArray = np.empty_like(pcSamples,dtype = float)
    
    phenoArray = np.empty_like(pcSamples,dtype = int)
    
def filter_pcs(iPath,lofString='hc_lof',f = phenoFile,pcPath = eigenvecPath):
    '''
    Filters the eigenvec file to keep only samples that are shared across all files
    '''
    samples = get_shared_samples(iPath,lofString,f, pcPath)
    print('samples loaded.')
    with open(pcPath,'rt') as i,open(iPath + lofString + '_pcs.txt','wt') as o:
        for line in i:
            sample = line.strip().split(' ')[0]
            if sample in samples:
                o.write(line)

def get_shared_samples(iPath,lofString = 'hc_lof',f = phenoFile,pcPath = eigenvecPath):
    '''
    Returns and saves the samples shared across all files
    '''
    sharedPath = dataPath +lofString + '_shared_samples.txt'
    if os.path.isfile(sharedPath):
        samples = np.loadtxt(sharedPath,dtype = str)
    else:
        print('importing all samples..')
        pcSamples = return_pc_samples(pcPath)
        lofSamples = return_lof_samples(iPath,lofString)
        phenoSamples = return_column(f =f,dtype = str)
        print('done')
        sampleLists = [pcSamples,lofSamples,phenoSamples]
        samples = set(sampleLists[0])
        for sampleList in sampleLists[1:]:
            samples = samples.intersection(sampleList)

        samples = np.array(list(samples))
        np.savetxt(sharedPath,samples,fmt ='%s')
    return samples

def return_lof_samples(iPath,lofString = 'hc_lof'):
    matrixPath =  iPath + lofString + "_gene_to_sample.tsv"
    with open(matrixPath,'rt') as i:
        samples =i.readline().strip().split('\t')[1:]
    return np.array(samples)

def return_pc_samples(pcPath = eigenvecPath):
 
    pcSamples = np.loadtxt(pcPath,dtype = str,usecols = [0])

    return pcSamples
    


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
