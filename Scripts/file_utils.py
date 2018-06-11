import os
import numpy as np
import gzip
from collections import defaultdict
import pickle

rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
dataPath = rootPath + 'Data/'
# REQUIRED FILES
rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

dataPath = rootPath + 'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = rootPath + 'tmp_scripts/'

def dd(tp):
    return defaultdict(tp)
def dd_str():
    return defaultdict(str)
    
def make_sure_path_exists(path):
    import errno
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise                

for path in [dataPath,bashPath]:
    make_sure_path_exists(path)

def return_header_variants(matrixPath):
    '''
    The plink command adds the alt to the name of the variant. Here i loop through the variants and just return the original variant name
    '''
    with open(matrixPath,'rt') as i:
        header = i.readline()
        header = header.strip().split(' ')[1:]
        headerVariants = ['_'.join(elem.split('_')[:-1]) for elem in header]
        
    return np.array(headerVariants,dtype = str)
                

def split_array_chunk(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


import shlex
from subprocess import call

def git_pull():
    call(shlex.split('git pull'))


def read_header(header = None,lofString = "hc_lof"):
    '''
    Reads the first line of the variants.gz file and returns the position of the Info score
    '''
    
    if header == None:
        with gzip.open(annotatedVariants,'rt') as i:
    
            header = i.readline()
            header = header.strip().split('\t')

            
    infoPos = [i for i,elem in enumerate(header) if 'INFO_' in elem]
    lofPos = [i for i,elem in enumerate(header) if  elem == lofString][0]
    avgPos = [i for i,elem in enumerate(header) if 'INFO' == elem][0]
    genePos =  [i for i,elem in enumerate(header) if  elem == "gene"][0]
    return infoPos,lofPos,avgPos,genePos

def return_column(pheno = 'FINNGENID',f = phenoFile,dtype = 'f8'):
    '''
    Given the phenotype file, it returns the column given a phenotype
    '''
    header = return_header(f =f )
    for i,elem in enumerate(header):
        if str(elem) == pheno:
            phenocol = i
    idcol = 0
    if f.split('.')[-1] == 'txt':
        i = open(f,'rt')
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


def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.



def get_variant_to_gene_dict(iPath,lofString = 'hc_lof'):
    '''
    Reads the plink snplist and returns a gene to variant dictionary 
    '''
    #get variant to gene mapping from full list of variants
    bFile = iPath + 'plink_files/'+lofString 

    v2g = dd(str)
    with open(dataPath + lofString + '_variants.txt','rt') as i:
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
    return g2v



def variant_is_dict(annVariants = annotatedVariants,iPath ='/home/pete/results/hc_lof/',lofString = "hc_lof" ):
    
    '''
    Read the annotated_variants and returns a dict[variant][batch] = INFO_SCORE for teh variants that are in the snplist.
    I can use this dictionary to retreieve the info score for the samples
    '''

    picklePath = iPath + '/misc/'
    make_sure_path_exists(picklePath)

    picklePath += lofString +'_vDict.p'
    print('loading/generating dict[variant][batch] = INFO_SCORE dict -->' + picklePath)
    try:
        vDict = pickle.load(open(picklePath,'rb'))
    except:
        snplist = iPath + '/plink_files/' +lofString + '.snplist'

        print('data missing, generating..')
        variants = np.loadtxt(snplist,dtype = str)   
        vDict = defaultdict(dd_str)
        with gzip.open(annVariants,'rt') as i:
            #read header
            header = i.readline().strip().split('\t')
            infoPos,lofPos,avgPos,genePos = read_header(header)
            #return position of batches 
            batches = header[infoPos[0]:infoPos[-1]+1]
            #return batches
            batches = [batch.split('INFO_')[1].split('_R1')[0] for batch in batches]
            pickle.dump(batches,open(dataPath + 'ourbatches.p','wb'))
        
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

        pickle.dump(vDict,open(picklePath,'wb'))

    return vDict



def sample_to_batch_dict(filePath = dataPath + 'sample_info.txt'):
    '''
    Given timo's file maps a sample to a batch. requires a conversion on the fly due to slightly different names between his batch names and ours. Need to pass our batches and use difflib
    '''

    picklePath = dataPath +  'sample_to_batch.p'

    try:
        print('returning sample to batch dict')
        s2b = pickle.load(open(picklePath,'rb'))
        
    except:
        print("data doesn't exist..generating..")
        import difflib                         
        ourbatches = pickle.load(open(dataPath + 'ourbatches.p','rb'))

                             
        s2b = dd(str)
        sampleData = np.loadtxt(filePath,dtype = str,delimiter=':',usecols=[0,-1])
        for entry in sampleData:
            timoBatch,sample = entry
            ourBatch = difflib.get_close_matches(timoBatch,ourbatches)[0]
            s2b[sample] = ourBatch

        pickle.dump(s2b,open(picklePath,'wb'))
    return s2b


