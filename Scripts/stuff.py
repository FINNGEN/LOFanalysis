



###############################################
#---CONVERTS LOF MATRIX TO INFO_SCORE MATIX---#
###############################################
def write_info_score_matrix(annotatedPath,oPath,lofString,batchPath = dataPath + 'sample_info.txt'):

    '''
    Goes through each line of the matrix(sample data) and updates the 1s to be the Info score for that sample's batch.
    process_line() is called for each line of the matrix
    '''
    
    oFile = oPath + lofString + '_info_score_matrix.txt'

    if os.path.isfile(oFile):
        print('info score matrix already generated')

    else:
        print('generating info score matrix..')
        matrixPath = oPath + '/plink_files/'+ lofString + matrixName

        #stuff required  
        vDict = variant_is_dict(annotatedPath,oPath,lofString)
        s2b = sample_to_batch_ditct(batchPath)
        headerVariants = return_header_variants(matrixPath)
    
        print('looping samples...')
        with open(matrixPath,'rt') as i,open(oFile,'wt') as o:
            next(i) #skip header
            for line in i:
                sample,data = process_line(line,s2b,headerVariants,vDict)
                o.write(sample + ' ' + ' '.join([str(elem) for elem in data]) + '\n')
    

def process_line(line,s2b,headerVariants,vDict):
    '''
    Given a line of the lof_matrix, it returns 0 if not lof and INFO_SCORE of the batch for that variant otherwise.
    It uses
    '''
    line = line.strip().split(' ')
    # get sample info
    sample = line[0]
    batch = s2b[sample]
    #read data from sample
    data = np.array(line[1:],dtype = str)
    data = np.isin(data,['1','2']).astype(float) # boolean if lof or not
    #now we have a 1 if there is lof and 0 elsewhere
    dataMask = np.where(data==1)[0] #index of variant with lof
    # now i create a mini array only for the positions with lof and i multiply the original array, masking it, with the info-score array
    infoArray = np.empty(len(dataMask),dtype = float)
    for i,elem in enumerate(infoArray):
        lofVariant = headerVariants[i]
        infoArray[i] = vDict[lofVariant][batch]
    data[dataMask] *= infoArray

    return sample,data
    



def variant_is_dict(annVariants = annotatedVariants,iPath ='/home/pete/results/hc_lof/',lofString = "hc_lof" ):
    
    '''
    Read the annotated_variants and returns a dict[variant][batch] = INFO_SCORE for teh variants that are in the snplist.
    I can use this dictionary to retreieve the info score for the samples
    '''

    picklePath = iPath + lofString + '_vDict.p'
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




def sample_to_batch_ditct(filePath = dataPath + 'sample_info.txt'):
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
