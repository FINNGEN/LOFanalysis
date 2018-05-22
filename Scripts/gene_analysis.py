import numpy as np
import os

rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
dataPath = rootPath + 'Data/'
annotatedVariants =  dataPath + 'annotated_variants.gz'
bashPath = rootPath + 'tmp_scripts/'
for path in [dataPath,bashPath]:
    make_sure_path_exists(path)



eigenvecPath = dataPath + '10pc.eigenvec'


def return_pc_samples(pcPath = eigenvecPath):

    
    pcSamples = np.loadtxt(pcPath,dtype = str,usecols = [0])

    return pcSamples
    
