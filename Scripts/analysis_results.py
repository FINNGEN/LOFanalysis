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
                for line in i:
                    o.write(gene + '\t' + line)
    
    cmd = 'cat ' +filePath + ' | sort -k8 > ' + iPath + lofString + '_gene_summary_ordered.txt'

    shPath = bashPath +  'filter_lof_matrix.sh'

    with open(shPath,'wt') as o:
        o.write(' #!/bin/bash\n')
        o.write(cmd)

    call(['chmod','+x',shPath])
    call(shPath,shell = True)
