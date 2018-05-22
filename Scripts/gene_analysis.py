import numpy as np
import os
import gzip

rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
dataPath = rootPath + 'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = rootPath + 'tmp_scripts/'

# REQUIRED FILES
phenoList = np.loadtxt(dataPath + 'pheno-list.txt',usecols = [0],dtype = str)
phenoFile = dataPath + 'FINNGEN_PHENOTYPES_DF1_V4_2018_03_27.txt.gz'
eigenvecPath = dataPath + '10pc.eigenvec'


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
