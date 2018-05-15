import numpy as np
import os
import gzip
from operator import itemgetter
import os
import shutil
from collections import defaultdict as dd
import shlex
from subprocess import Popen, PIPE,call

plink = shutil.which('plink')

currPath = os.getcwd() + '/'
dataPath =  'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = 'tmp_scripts/'


def get_variant_to_gene_dict():

    v2g = dd(str)
    with open(dataPath + 'lof_variants.txt','rt') as i:
        for line in i:
            variant,gene = line.strip().split('\t')
            v2g[variant] = gene

    return v2g



def generate_matrix(iFile,oFile):
    """
    Returns variant x sample matrix with 1s where variant is present
    """
    cmd = 'plink -bfile '+ iFile +' --recode A --out ' + oFile
    call(shlex.split(cmd))
    #remove unncessary columns
    cmd = "cat " +oFile + " |cut -d ' ' -f-2,7- > lof_matrix.tsv "
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

