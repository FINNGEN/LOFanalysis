import numpy as np
import os
import gzip
from collections import defaultdict as dd
import sys

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


#####################################
#--FIX FILES TO ORDER SAMPLE DATA---#
#####################################

def transpose_matrix(iPath,lofString= 'hc_lof'):
    import csv
    import itertools


    # Read all Tab-delimited rows from stdin.
    oMatrix = iPath + lofString + "_gene_to_filtered_sample.tsv"

    tsv_reader = csv.reader(open(oMatrix,'rt'), delimiter='\t')
    all_data = list(tsv_reader)
    print('data imported..')
        
    # Transpose it.
    all_data = list(itertools.zip_longest(*all_data, fillvalue=''))
    print('data transposed...')
    # Write it back out.
    tsv_writer = csv.writer(open(oMatrix +'.tmp','wt'), delimiter='\t')
    for row in all_data:
        tsv_writer.writerow(row) 

def reorder_lof_matrix(iPath,lofString = 'hc_lof'):
    '''
    Shuffles the column of the gene_lof_matrix so that the ordering of samples is the same as in the eigenvec file
    '''
    

    iMatrix = iPath + lofString + "_gene_to_sample.tsv"
    oMatrix = iPath + lofString + "_gene_to_filtered_sample.tsv"

    if os.path.isfile(oMatrix):
        print('matrix already reordered')
    else:
        samples = get_shared_samples(iPath,lofString)
        lofSamples = return_lof_samples(iPath,lofString)
        with open(iMatrix,'rt') as i,open(oMatrix ,'wt') as o:
            for line in i:
                # i create a dict that stores the info for each samples so i can then rearrange them
                geneDict = dd()
                line = line.strip().split('\t')
                gene = line[0]
                sys.stdout.write('processing gene %s \r'%(gene)),
                sys.stdout.flush()
                data = line[1:]
                assert len(data) == len(lofSamples)
                # store sample Data
                for j,sample in enumerate(lofSamples):
                    geneDict[sample] = data[j]
                # write new line only keeping shared sample data
                newLine = gene +'\t'
                newLine += '\t'.join([geneDict[tmpSample] for tmpSample in samples])
                o.write(newLine + '\n')


        
def return_pc_samples(iPath,lofString = 'hc_lof'):

    pcPath = iPath + lofString + '_pcs.txt'
    pcSamples = np.loadtxt(pcPath,dtype = str,usecols = [0])

    return pcSamples
    
    
def filter_pcs(iPath,lofString='hc_lof',f = phenoFile,pcPath = eigenvecPath):
    '''
    Filters the eigenvec file to keep only samples that are shared across all files
    '''

    pcFile = iPath + lofString + '_pcs.txt'
    print('filtering pcs to shared samples --> ' + pcFile)
    if os.path.isfile(pcFile):
        print('principal components already filtered')
    else:
        samples = get_shared_samples(iPath,lofString,f, pcPath)
        print('samples loaded.')
        with open(pcPath,'rt') as i,open(pcFile,'wt') as o:
            for line in i:
                sample = line.strip().split(' ')[0]
                if sample in samples:
                    o.write(line)




###############################
#----RETURN SHARED SAMPLES----#
###############################



def get_shared_samples(iPath,lofString = 'hc_lof',f = phenoFile,pcPath = eigenvecPath):
    '''
    Returns and saves the samples shared across all files and it returns them in the order of the pc file.
    After the first reordering, the function will always load the same order of samples
    '''
    sharedPath = iPath +lofString + '_shared_samples.txt'
    if os.path.isfile(sharedPath):
        finalSamples = np.loadtxt(sharedPath,dtype = str)
    else:
        print('importing all samples..')
       

        def return_pc_original_samples(pcPath = eigenvecPath):
            return np.loadtxt(pcPath,dtype = str,usecols = [0])
       
        pcSamples = return_pc_original_samples(pcPath)
        lofSamples = return_lof_samples(iPath,lofString)
        phenoSamples = return_column(f =f,dtype = str)
        print('done')
        sampleLists = [pcSamples,lofSamples,phenoSamples]
        samples = set(sampleLists[0])
        for sampleList in sampleLists[1:]:
            samples = samples.intersection(sampleList)

        samples = np.array(list(samples))
        print('set of samples calculated')
        finalSamples = [s for s in pcSamples if s in samples]
        np.savetxt(sharedPath,finalSamples,fmt ='%s')
    return finalSamples


def return_lof_samples(iPath,lofString = 'hc_lof'):
    matrixPath =  iPath + lofString + "_gene_to_sample.tsv"
    with open(matrixPath,'rt') as i:
        samples =i.readline().strip().split('\t')[1:]
    return np.array(samples)

##############################
#--PARSE THE PHENOTYPE FILE--#
##############################

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



if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description="Analyze lof data")
    parser.add_argument("--phenoFile", type= str,
                        help="path to phenotype file",required = False,default = phenoFile)
    parser.add_argument("--pcPath", type= str,
                        help="path to eigenvec file",required = False,default = eigenvecPath)



    subparsers = parser.add_subparsers(help='help for subcommand',dest ="command")

    # create the parser for the generate_variants command
    parser_fix_samples = subparsers.add_parser('fix-samples', help='fix files in order to match shared samples')
    parser_fix_samples.add_argument("--lof", type= str,
                        help="type of lof filter",required = True )
    parser_fix_samples.add_argument("--oPath", type= str,help="Path to folder where to output",default = ".")
    # create the parser for the "command_2" command
    
    
    args = parser.parse_args()

    if args.command == "fix-samples":
        oPath = (args.oPath + '/' + args.lof +'/').replace('//','/')
        filter_pcs(oPath,args.lof,args.phenoFile,args.pcPath)

    
    
