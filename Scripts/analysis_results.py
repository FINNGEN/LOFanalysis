import numpy as np
import os
from operator import itemgetter
import pickle
import shlex
import sys
from subprocess import Popen, PIPE,call
from file_utils import make_sure_path_exists,gzip,get_filepaths
from file_utils import rootPath,dataPath,annotatedVariants,bashPath
import matplotlib as mpl
from matplotlib import pyplot as plt
from verkko.plots import matplotlibHelperFunction as HF
import pylab

miscPath = rootPath+ '/misc/'


def qq_data(resPath ,lofString = 'hc_lof'):

    resPath += '/' + lofString +'/results/' 
    files = get_filepaths(resPath)
    for f in files[:10]:
        with gzip.open(f,'rt') as i:
            count = 0
            while count <10:
                for line in i:
                    pval = line.split('\t')[1]
                    try:
                        pval = float(pval)
                    except:
                        pass
                    print(pval)
                    count +=1
    
    

def qq_plot(resPath ,lofString = 'hc_lof',figPath = '/home/pete/results/'):

    figPath += lofString + '/qq_plot.pdf'
    pylab.ioff()
    fig = HF.setFigure()
    gs = mpl.gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0,0])
    print('ax created...')
#    plt.legend(scatterpoints=1, frameon=False,labelspacing=1, loc='lower left');

    print('saving...')
    fig.savefig(figPath)
    plt.close(fig)

def scatter_file(lofString = 'hc_lof',infoFilter = 0.9):

    resultsPath = miscPath + '/tmp.txt'
    spaResults = '/home/pete/Data/SPA_data/SPA_results'
    #FISHER RESULTS
    data = np.loadtxt(resultsPath,dtype = str, delimiter = '\t',usecols=[0,1])
    genes = data[:,0]
    phenotypes = data[:,1]
    fisherResults = np.loadtxt(resultsPath,dtype = float, delimiter = '\t',usecols=[-2,-1])
    # CREATE CORRELATION
    pvals = np.empty((len(genes),2),dtype = float)
    oddsratio = np.empty_like(pvals)
    for j,gene in enumerate(genes):
        pheno = phenotypes[j]
        print(gene,pheno)
        phenoResults = spaResults+'/'+ pheno + '-' + lofString + '_gene_to_sample_'+str(infoFilter) + '.txt.gz'
        with gzip.open(phenoResults,'rt') as i:
            for line in i:
                line = line.split('\t')
                if line[0] == gene:
                    pval = float(line[1])
                    odd = float(line[-2])
                    print(fisherResults[j])
                    pvals[j] = [fisherResults[j,1],pval]
                    oddsratio[j] = [fisherResults[j,0],odd]

    np.savetxt(miscPath + 'pvals.txt',pvals,fmt = '%E')
    np.savetxt(miscPath + 'odds.txt',oddsratio,fmt = '%f')
    return None

                
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
                    pval = float(line[-1])
                    if (odds > 1) & (pval < 0.00000001):
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

    genes = keep_single_pheno(iPath,lofString)
    print(genes)
    with open(oFile,'rt') as i :
        resLines = []
        next(i)
        for line in i:
            line = line.split('\t')
            gene = line[0]
            print(gene)
            odds = line[-2]
            line[-2] = str(round(float(odds),5))
            pval = line[-1]
            line[-1] = format_e(float(pval))
            line = '\t'.join(line)
            if gene in genes:
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
    print(oFile)
    genes = np.loadtxt(oFile,dtype = str,usecols = [0])
    from collections import Counter
    c = Counter(genes)
    return [gene for gene in c if c[gene] ==1]
