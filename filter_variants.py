import numpy as np
import os
import gzip
from operator import itemgetter
import os
import shutil
from collections import defaultdict 

import shlex
from subprocess import Popen, PIPE,call

plink = shutil.which('plink')

currPath = os.getcwd() + '/'
dataPath =  'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = 'tmp_scripts/'

def dd():
    return defaultdict(int)

def sample_to_batch_ditct(filePath):
    '''
    Given timo's file maps a sample to a batch
    '''
    s2b = dd()
    with open(filePath,'rt') as i:
        for line in i:
            line = line.strip().split(':')
            sample = line[-1]
            batch = line[0]
            s2b[sample] = batch
    return s2b

def variant_is_dict(snplist ='/home/pete/lof_data/filtered_lof.snplist' ):
    
    '''
    Read the annotated_varaints and returns a dict[variant][batch] = INFO_SCORE
    '''

    variants = np.loadtxt(snplist,dtype = str)
    vDict = defaultdict(dd)
    
    with gzip.open(annotatedVariants,'rt') as i:
        header = i.readline()
        infoPos,lofPos,avgPos,genePos = read_header(header.strip().split('\t'))
        
        
def return_gene_columns(gene,filePath,g2v):
    """
    Loops through the header of the matrix file and returns the columns where variants belong to the gene
    """
    variantList = g2v[gene]
    with open(filePath,'rt') as i:
        header = i.readline()
        header = header.strip().split(' ')
    geneColumns = [i for i,elem in enumerate(header) if '_'.join(elem.split('_')[:-1]) in variantList]

    #import sample data keeping columns of gene
    vData = np.loadtxt(filePath,dtype = str,usecolds = geneColumns,skiprows = 1)
    #convert NA to 0
    vData[vData =='NA'] = 0
    #convert to int
    vData = vData.astype(int)
    #sum across variants and check if >1
    gData = (np.sum(sampleData,axis = 1) >0).astype(int)
    return gData


def get_variant_to_gene_dict(bFile):

    #get variant to gene mapping from full list of variants
    v2g = dd(str)
    with open(dataPath + 'lof_variants.txt','rt') as i:
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
    return v2g,g2v




#######################
#--GENERATING MATRIX--#
#######################
def generate_matrix(iFile,oFile):
    """
    Returns variant x sample matrix with 1s where variant is present
    """
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


            
def plink_filter(filePath,snpslist,oPath,geno = 0.9):
    """
    Filter full data for only varianst we need
    """
    
    cmd = 'plink -bfile ' + filePath + ' --geno ' + str(geno) + ' --extract ' + snpslist + ' --make-bed -out ' + oPath
    call(shlex.split(cmd))
    cmd = 'plink -bfile ' + oPath + ' --write-snplist --out ' + oPath
    call(shlex.split(cmd))
    
def create_info_file():

    '''
    Creates a lof_variants.txt with variants that carry lof along with their genes
    '''
    
    with gzip.open(annotatedVariants,'rt') as i,open(dataPath + 'lof_variants.txt','wt') as o:
        infoPos,lofPos,avgPos,genePos = read_header(i.readline().strip().split('\t'))

        for line in i:
            line = line.strip().split('\t')
            variant = line[0]
            lof = line[lofPos]
            gene = line[genePos]
            if (lof == "true"):
                o.write(variant.replace(':','_') + '\t' + gene + '\n')

    #write snplist for plink
    shPath = bashPath +  'snplist.sh'
    cmd = "cat Data/lof_variants.txt | cut -f1 >> Data/lof.snplist"
    with open(shPath,'wt') as o:
        o.write(' #!/bin/bash\n')
        o.write(cmd)

    call(['chmod','+x',shPath])
    call(shPath,shell = True)


def read_header(header = None):
    '''
    Reads the first line of the variants.gz file and returns the position of the Info score
    '''
    
    if header == None:
        with gzip.open(annotatedVariants,'rt') as i:
    
            header = i.readline()
            header = header.strip().split('\t')

            
    infoPos = [i for i,elem in enumerate(header) if 'INFO_' in elem]
    lofPos = [i for i,elem in enumerate(header) if  elem == "hc_lof"][0]
    avgPos = [i for i,elem in enumerate(header) if 'INFO' == elem][0]
    genePos =  [i for i,elem in enumerate(header) if  elem == "gene"][0]
    return infoPos,lofPos,avgPos,genePos

