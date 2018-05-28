import numpy as np
import os
from operator import itemgetter
import pickle
import shlex
import sys
from subprocess import Popen, PIPE,call
from file_utils import make_sure_path_exists,return_header_variants,split_array_chunk,read_header,get_filepaths


rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'


def write_final_file(iPath,lofString = 'hc_lof'):

    oPath = iPath + '/gene_fits/'
    fileList = [p for p in get_filepaths(oPath) if p.endswith('_gene_results.txt')]
    filePath = iPath + lofString + '_gene_summary.txt'
    with open(filePath,'wt') as o:
        for f in fileList:
            gene = f.split('most_severe')[-1].split('_')[1]
            print(gene,filePath)
            with open(f,'rt') as i:
                for line in i:
                    o.write(gene + '\t' + line)
    
