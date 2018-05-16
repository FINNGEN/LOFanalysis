import numpy as np
import os
import gzip
from operator import itemgetter
import shutil
from collections import defaultdict 
import pickle
import shlex
from subprocess import Popen, PIPE,call
from file_utils import make_sure_path_exists
plink = shutil.which('plink')

currPath = os.getcwd() + '/'
dataPath =  'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = 'tmp_scripts/'
for path in [dataPath,bashPath]:
    make_sure_path_exists(path)


def dd(tp):
    return defaultdict(tp)

lofName = "filtered_lof"
matrixName = "lofvariantmatrix.tsv"


###############################
#--MERGE VARIANTS INTO GENES--#
###############################

def write_new_matrix(iPath):

    g2v = get_variant_to_gene_dict(iPath)
    samples =  np.loadtxt(iPath + matrixName,dtype = str,usecols =[0])
    with open(iPath + "gene_to_sample_lof.txt",'wt') as f:
        f.write("\t".join(samples) + '\n')
        for gene in g2v:
            gData = return_gene_columns(gene,iPath,g2v).astype(str)
            gArray = np.concatenate((np.array([gene]),gData))
            assert gArray.shape == samples.shape
            f.write("\t".join(gArray) + '\n')
    

def return_gene_columns(gene,iPath,g2v):
    """
    Loops through the header of the matrix file and returns the columns where variants belong to the gene
    """
    geneVariants = g2v[gene]
    filePath = iPath + matrixName
    #TEST
    headerVariants = return_header_variants(filePath)
    geneColumns = [i+1 for i,elem in enumerate(headerVariants) if elem in geneVariants]
    print(geneColumns)

    #import sample data keeping columns of gene
    vData = np.loadtxt(filePath,dtype = str,usecols = geneColumns,skiprows = 1)
    #convert NA to 0
    vData[vData =='NA'] = 0
    #convert to int    
    vData = vData.astype(int)        
    if len(geneColumns) > 1:
        #sum across variants and check if >1
        vData = np.sum(vData,axis = 1)
    
    return (vData>0).astype(int)


def get_variant_to_gene_dict(iPath):
    '''
    Reads the plink snplist and returns a gene to variant dictionary 
    '''
    #get variant to gene mapping from full list of variants
    bFile = iPath + lofName

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
    return g2v

def return_header_variants(iPath):
    '''
    The plink command adds the alt to the name of the variant. Here i loop through the variants and just return the original variant name
    '''
    with open(iPath + matrixName,'rt') as i:
        header = i.readline()
        header = header.strip().split(' ')[1:]
        headerVariants = ['_'.join(elem.split('_')[:-1]) for elem in header]
        
    return np.array(headerVariants,dtype = str)
                
                      



#######################
#--GENERATING MATRIX--#
#######################
def generate_matrix(iPath):
    """
    Returns variant x sample matrix with 1s where variant is present
    """
    oFile = iPath + matrixName
    iFile = iPath + lofName
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


            
def plink_filter(filePath,oPath,geno = 0.9):
    """
    Filter full data for only varianst we need
    """
    snpslist = dataPath + "lof.snplist"
    make_sure_path_exists(oPath)
    cmd = 'plink -bfile ' + filePath + ' --geno ' + str(geno) + ' --extract ' + snpslist + ' --make-bed -out ' + oPath + lofName
    call(shlex.split(cmd))
    cmd = 'plink -bfile ' + oPath + lofName +  ' --write-snplist --out ' + oPath + lofName
    call(shlex.split(cmd))

    generate_matrix(oPath)
    
def create_info_file(annotatedFile):
    '''
    Creates a lof_variants.txt with variants that carry lof along with their genes
    '''
    lofPath =dataPath + 'lof_variants.txt'
    if os.path.isfile(lofPath):
        print('variants already filtered')
        return 

    else:
        with gzip.open(annotatedFile,'rt') as i,open(lofPath,'wt') as o:
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



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Generation of gene to samples LOF matrix")

    
    parser.add_argument("--annotatedFile", type= str,
                        help="path to annotatedFile",required = False,default =annotatedVariants )
    parser.add_argument("--plinkFile",type = str,help ="path to plink files")
    parser.add_argument("--oPath",type = str,help ="file to write results")
    parser.add_argument("--geno",type = float,required = False,default = 0.9,help ="genotype call rate")
    args = parser.parse_args()
   
    create_info_file(args.annotatedFile)
    plink_filter(args.plinkFile,args.oPath,geno = args.geno)
    
