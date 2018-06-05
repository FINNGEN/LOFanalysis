def multiproc_logit_gene(iPath,lofString='hc_lof',f = phenoFile,proc = cpus,test = True,infoFilter = 0.9):

    geneList = get_info_score_gene_list(iPath,lofString,infoFilter)
    print(len(geneList),' genes')
    geneList= geneList if test is False else geneList[:proc]
        

    params  = list(product([iPath],geneList,[lofString],[f],[test]))
    pool = multiprocessing.Pool(proc)
    pool.map(logit_wrapper_gene,params)
    pool.close()

def logit_wrapper_gene(args):
    logistic_gene(*args)
    
def logistic_gene(iPath,gene,lofString = 'hc_lof',f = phenoFile,test = True,infoFilter = 0.9):
    '''
    Function that is ultimately passed to the multiprocessing pool. It loops through all genes given a phenotype. With test  it only works with a small chunk of genes
    '''
    oPath = iPath + '/fits/'
    make_sure_path_exists(oPath)
    oFile = oPath + lofString + '_' + gene + '_' + str(infoFilter) + '_gene_results.txt'
 
    print(gene)
    lofData = get_lof_data(iPath,gene,lofString)

 
    pcPath = iPath + lofString + '_pcs.txt'
    pcData = np.loadtxt(pcPath,dtype = float,usecols = range(1,11))
    print('pcData loaded.')
    
    with open(oFile,'wt') as o:
        o.write('\t'.join(shlex.split('gene lof_cases lof_controls no_lof_cases no_lof_controls logit_coeff_gene logit_pval_gene logit_coeff_pc1 logit_pval_pc1 logit_coeff_pc2 logit_pval_pc2 fischer_oddsratio fischer_pval ')) + '\n')

        if test is True:
            pList = phenoList[:20]
        print(len(pList))

        for i,pheno in enumerate(pList):
            #print(pheno,i,gene)
            o.write(pheno + '\t')
            phenoDataPath = iPath + '/pheno_data/'
            make_sure_path_exists(phenoDataPath)
            phenoSave = phenoDataPath + lofString + '_' + pheno  + '_phenodata.txt'
            try:
                phenoData = np.loadtxt(phenoSave,dtype = int)
            except:
                phenoData = get_pheno_data(iPath,pheno,f,lofString)
                np.savetxt(phenoSave,phenoData,fmt = '%i')

            logit_results,f_results,table = logistic_regression(iPath,lofString,pcData,phenoData,lofData,f)
            # write counts of lof/no_lof
            countString =  '\t'.join([str(elem) for elem in table.flatten()])
            o.write(countString + '\t')
            
            try:
                 params = logit_results.params
                 pvalues = logit_results.pvalues
                 #add columns
                 res = np.column_stack((params,pvalues))
                 # flatten so first two elemts are from lof, next 2 pc1 etc.
                 resArray = res.flatten()[:6]
            except:
                resArray = ['NA']*6
                
            # write logit_results            
            oString =  '\t'.join([str(elem) for elem in resArray])
            o.write( oString + '\t')
            o.write( '\t'.join([str(elem) for elem in f_results]))
            o.write('\n')
    return None



   
def multiproc_write_pheno(iPath,lofString='hc_lof',f = phenoFile,proc = cpus,test = True):

    pList = phenoList if test is False else phenoList[:proc]
    print(len(pList))
    params  = list(product([iPath],pList,[lofString],[f],[test]))
    pool = multiprocessing.Pool(proc)
    pool.map(pheno_wrapper,params)
    pool.close()

def pheno_wrapper(args):
    return_pheno_data(*args)

##############################
#--PARSE THE PHENOTYPE FILE--#
##############################

def return_pheno_data(iPath,pheno,lofString = 'hc_lof',f = phenoFile,test = True):
    oPath = iPath + '/pheno_data/'
    make_sure_path_exists(oPath)
    oFile = oPath + lofString + '_' + pheno  + '_phenodata.txt'
    print(pheno)
    phenoData = get_pheno_data(iPath,pheno,f,lofString)
    np.savetxt(oFile,phenoData,fmt = '%i')
    return None
    


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
    

def sample_to_batch_ditc(filePath = dataPath + 'sample_info.txt'):
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




