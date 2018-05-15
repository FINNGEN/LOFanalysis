import numpy as np
import os
import gzip
from operator import itemgetter

import pickle
currPath = os.getcwd() + '/'
dataPath =  'Data/'
filePath =  dataPath + 'annotated_variants.gz'

from collections import defaultdict as dd

import shlex
from subprocess import Popen, PIPE



def generate_vars_file():

    '''
    Returns the list of positions for tabix from the unzipped filtered_variants.file
    '''
    #keeps variants and minimum score
    cmd = """gunzip -c filtered_variants.gz | awk 'BEGIN {OFS = ","} {print $1,$3}' """
    # keeps only the variants and remove header
    cmd = """gunzip -c filtered_variants.gz | head | cut -f1 |sed 1d"""



def create_info_file():

    '''
    Creates a filtered_variants_score.txt file which only contains the variants with minimum IS across batches >= filt
    '''
    
    with gzip.open(filePath,'rt') as i,open(dataPath + 'lof_variants.txt','wt') as o:
        infoPos,lofPos,avgPos,genePos = read_header(i.readline().strip().split('\t'))

        for line in i:
            line = line.strip().split('\t')
            variant = line[0]
            info = min([float(elem) for elem in itemgetter(*infoPos)(line)])
            avg = float(line[avgPos])
            lof = line[lofPos]
            gene = line[genePos]
            if (lof == "true"):
                o.write(variant.replace(':','_') + '\t' + gene + '\n')

def read_header(header = None):
    '''
    Reads the first line of the variants.gz file and returns the position of the Info score
    '''
    
    if header == None:
        with gzip.open(filePath,'rt') as i:
    
            header = i.readline()
            header = header.strip().split('\t')

            
    infoPos = [i for i,elem in enumerate(header) if 'INFO_' in elem]
    lofPos = [i for i,elem in enumerate(header) if  elem == "hc_lof"][0]
    avgPos = [i for i,elem in enumerate(header) if 'INFO' == elem][0]
    genePos =  [i for i,elem in enumerate(header) if  elem == "gene"][0]
    return infoPos,lofPos,avgPos,genePos

