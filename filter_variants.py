import numpy as np
import os
import gzip
from operator import itemgetter
import shutil
from collections import defaultdict 
import pickle
import shlex
from subprocess import Popen, PIPE,call

plink = shutil.which('plink')

currPath = os.getcwd() + '/'
dataPath =  'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = 'tmp_scripts/'

timo_to_pietro_dict = {''}
def dd_str():
    return defaultdict(str)
def dd(tp):
    return defaultdict(tp)


###############################
#--MERGE VARIANTS INTO GENES--#
###############################

def write_new_matrix(g2v,filePath,oFile):

    samples =  np.loadtxt(filePath,dtype = str,usecols =[0])
    fileHandle = open(oFile,'a')
    np.savetxt(fileHandle,samples,fmt = '%s')
    for gene in g2v:
        gData = return_gene_columns(gene,filePath,g2v).astype(str)
        gArray = np.concatenate((np.array([gene]),gData))
        assert gArray.shape == samples.shape
        np.savetxt(fileHandle,gArray,fmt = '%s',newline =" ")
        break
    fileHandle.close()

def return_gene_columns(gene,filePath,g2v):
    """
    Loops through the header of the matrix file and returns the columns where variants belong to the gene
    """
    geneVariants = g2v[gene]
  
    #TEST
    headerVariants = return_header_variants(filePath)
    geneColumns = [i+1 for i,elem in enumerate(headerVariants) if elem in geneVariants]
    print(geneColumns)

    #import sample data keeping columns of gene
    vData = np.loadtxt(filePath,dtype = str,usecols = geneColumns,skiprows = 1)
    
    if len(geneColumns) > 1:
        
        #convert NA to 0
        vData[vData =='NA'] = 0
        #convert to int
        vData = vData.astype(int)
        #sum across variants and check if >1
        return (np.sum(vData,axis = 1) >0).astype(int)
    else:
        return vData


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



####################################
#--GET INFO SCORE OF LOF VARIANTS--#
####################################

def write_info_score_matrix(annotatedPath,snplist,batchPath,matrixPath,oPath):
    s2b = sample_to_batch_ditct(batchPath)
    vDict = variant_is_dict(snplist)
    headerVariants = return_header_variants(matrixPath)
    with open(matrixPath,'rt') as i,open(oPath,'wt') as o:
        next(i) #skip header
        for line in i:
            sample,data = process_line(line,s2b,headerVariants,vDict)
            o.write(sample + '\t' + '\t'.join([str(elem) for elem in data]) + '\n')
    
              
            

    

def process_line(line,s2b,headerVariants,vDict):
    '''
    Given a line of the lof_matrix, it returns 0 if not lof and INFO_SCORE of the batch for that variant otherwise
    '''
    line = line.strip().split(' ')
    sample = line[0]
    batch = s2b[sample]
    data = np.array(line[1:],dtype = str)
    data = np.isin(data,['1','2']).astype(float)
    #now we have a 1 if there is lof and 0 elsewhere
    dataMask = np.where(data==1)[0] #index of variant with lof
    # now i create a mini array that will multiply the 1s
    infoArray = np.empty(len(dataMask),dtype = float)
    for i,elem in enumerate(infoArray):
        lofVariant = headerVariants[i]
        infoArray[i] = vDict[lofVariant][batch]
    data[dataMask] *= infoArray

    return sample,data
    


def sample_to_batch_ditct(filePath):
    '''
    Given timo's file maps a sample to a batch. requires a conversion on the fly due to slightly different names between his batch names and ours. Need to pass our batches and use difflib
    '''
    s2b = dd(str)
    with open(filePath,'rt') as i:
        for line in i:
            line = line.strip().split(':')
            sample = line[-1]
            batch = line[0]
            s2b[sample] = batch
    return s2b


#---VARIANT/SAMPLE/INFO_SCORE DICT--#
def variant_is_dict(snplist ='/home/pete/lof_data/filtered_lof.snplist' ):
    
    '''
    Read the annotated_variants and returns a dict[variant][batch] = INFO_SCORE for teh variants that are in the snplist
    '''

    try:
        print('pickling..')
        vDict = pickle.load(open(dataPath + 'vDict.p','rb'))
    except:
        print('data missing, generating..')
        variants = np.loadtxt(snplist,dtype = str)   
        vDict = defaultdict(dd_str)
        with gzip.open(annotatedVariants,'rt') as i:
            header = i.readline().strip().split('\t')
            infoPos,lofPos,avgPos,genePos = read_header(header)
            batches = header[infoPos[0]:infoPos[-1]+1]
            batches = [batch.split('INFO_')[1].split('_R1')[0] for batch in batches]
            startPos = infoPos[0]
            rangebatches = np.arange(len(batches))
            assert len(batches) == len(infoPos)

            #loop variants
            for line in i:
                line = line.strip().split('\t')
                variant = line[0].replace(':','_')
                if variant in variants:
                    for b in rangebatches:
                        batch = batches[b]
                        vDict[variant][batch] = line[startPos + b]

        pickle.dump(vDict,open(dataPath + 'vDict.p','wb'))

    return vDict




def return_header_variants(filePath):
    '''
    The plink command adds the alt to the name of the variant. Here i loop through the variants and just return the original variant name
    '''
    with open(filePath,'rt') as i:
        header = i.readline()
        header = header.strip().split(' ')[1:]
        headerVariants = ['_'.join(elem.split('_')[:-1]) for elem in header]
        
    return np.array(headerVariants,dtype = str)
                
                      





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

