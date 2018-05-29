import numpy as np
import os
from operator import itemgetter
import pickle
import shlex
import sys
from subprocess import Popen, PIPE,call
from file_utils import make_sure_path_exists,return_header_variants,split_array_chunk,read_header,get_filepaths


rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
bashPath = rootPath + 'tmp_scripts/'


def write_final_file(iPath,lofString = 'hc_lof'):

    oPath = iPath + '/gene_fits/'
    fileList = [p for p in get_filepaths(oPath) if p.endswith('_gene_results.txt')]
    filePath = iPath + lofString + '_gene_summary.txt'
    with open(filePath,'wt') as o:
        for f in fileList:
            gene = f.split(lofString)[-1].split('_')[1]
            print(gene)
            with open(f,'rt') as i:
                next(i)
                for line in i:
                    line = line.split('\t')
                    odds = float(line[-2])
                    pval = float(line[1])
                    if (odds > 1) & (pval < 0.0000000001):
                        line = '\t'.join(line)
                        o.write(gene + '\t' + line )
                   # line[-2] = str(round(float(odds),5))
                   # pval = line[-1]
                   # line[-1] = format_e(float(pval))
                   # line = '\t'.join(line)


    oFile =  iPath + lofString + '_gene_summary_ordered.txt'
    cmd = '  sort -g -k8 '+filePath+ ' > ' + oFile + ' && head -n11 ' + oFile + ' > '+iPath+'temp.txt'

    shPath = bashPath +  'gene_results.sh'

    with open(shPath,'wt') as o:
        o.write(' #!/bin/bash\n')
        o.write(cmd)

    call(['chmod','+x',shPath])
    call(shPath,shell = True)
    
    with open(iPath + 'temp.txt','rt') as i :
        resLines = []
        next(i)
        for line in i:
            line = line.split('\t')
            odds = line[-2]
            line[-2] = str(round(float(odds),5))
            pval = line[-1]
            line[-1] = format_e(float(pval))
            line = '\t'.join(line)
            resLines.append(line)
    with open(iPath + 'temp.txt','wt') as o :
        o.write("gene\tpheno\tlof_cases\tlof_controls\tno_lof_cases\tno_lof_controls\todds ratio\tpval\n")
        for line in resLines:
            o.write(line + '\n')
def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]


def keep_single_pheno(iPath,lofString = 'hc_lof'):

    oFile =  iPath + lofString + '_gene_summary_ordered.txt'
    genes = np.loadtxt(oFile,dtype = str,usecols = [0])
    return genes
