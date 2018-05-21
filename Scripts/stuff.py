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




