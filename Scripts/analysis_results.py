import numpy as np
import os
from operator import itemgetter
import pickle
import shlex
import sys
from subprocess import Popen, PIPE,call
from file_utils import make_sure_path_exists,gzip,get_filepaths,compute_qq,gc_value
from file_utils import rootPath,dataPath,annotatedVariants,bashPath
import matplotlib as mpl
from matplotlib import pyplot as plt
from verkko.plots import matplotlibHelperFunction as HF
import pylab
import pandas as pd

miscPath = rootPath+ '/misc/'
figPath = rootPath + '/Figs/'


def best_hits(resPath ,lofString = 'hc_lof',exp = 6):
    qqPath = resPath + '/' + lofString+ '/' +lofString + '_qq_data.txt'
    oPath =  resPath + '/' + lofString+ '/' +lofString + '_best_hits.txt'
    resPath += '/' + lofString +'/results/'
    print('fetching data from ' + resPath)
    files = get_filepaths(resPath )
    for f in files:
        pheno = f.split('/')[-1].split('-')[0]
        print(pheno)
        with gzip.open(f,'rt') as i,open(oPath,'wt') as o:
            for line in i:
                line = line.split('\t')
                pval = line[1]
                gene = line[0]

                p = np.float128(pval)
                pExp = -np.log10(p)
                print(pExp)
                if pExp > exp:
                    oString = '\t'.join([pheno,gene,pval])
                    o.write(oString + '\n')

def qq_data(resPath ,lofString = 'hc_lof'):

    qqPath = resPath + '/' + lofString+ '/' +lofString + '_qq_data.txt'
    resPath += '/' + lofString +'/results/'
    print('fetching data from ' + resPath)
    files = get_filepaths(resPath )
    res = []
    for f in files:
        print(f)
        with gzip.open(f,'rt') as i:
            for line in i:
                pval = line.split('\t')[1]
                try:
                    pval = np.float128(pval)
                    res.append(pval)
                except:
                    pass
                
    res = np.array(res,dtype = np.float128)
    np.savetxt(qqPath,res,fmt = '%E')
    return res

def return_gc(qqPAth,lofString):
    qqData = pd.read_csv(qqPath,dtype =float,header = None).values.flatten()
    qqData.sort()
    qqData = qqData[qqData >0]
    return gc_value(qqData,0.5)

def qq_plot(qqPath,fPath = figPath,lofString = 'hc_lof',dpi = 300,nBins = 1000):

   
    fPath += lofString + '_qq_plot.pdf'
    qqData = pd.read_csv(qqPath,dtype =float,header = None).values.flatten()
    qqData.sort()
    qqData = qqData[qqData >0]
    print(gc_value(qqData,0.5))
    qqData = np.log10(qqData)*-1
    pylab.ioff()
    fig = HF.setFigure()
    gs = mpl.gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0,0])

    print('ax created...')
    #    entries = np.power(10,max(qqData))
    plotData = np.array(compute_qq(qqData,nBins))
    xData = plotData[:,0]
    yData = plotData[:,1]

    ax.plot(xData,yData,'bo')
    ax.set_xlim([0,xData.max()])
    ax.set_ylim([0,yData.max()])
    ax.set_aspect('equal')
    print('saving...')
    fig.savefig(fPath,dpi = dpi)
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
