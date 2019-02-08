from file_utils import *

r2_blacklist = ['CD2_BENIGN_LEIOMYOMA_UTERI', 'GEST_DIABETES', 'PRE_OR_ECLAMPSIA', 'CD2_BENIGN_OVARY', 'E4_OVARDYS', 'ENDOMETRIOSIS_']

def merge_gender(res_path):

    all_files = get_filepaths(res_path)
    print(len(all_files))
    res_files = []
    for f in all_files:
        if 'MALE.SAIGE.' in f:
            res_files.append(f)

    print(len(res_files))
    for f in all_files:
        if f not in res_files:
            #check if male specific phenotype exists
            male_test = f.replace('.SAIGE','.MALE.SAIGE')
            female_test = f.replace('.SAIGE','.FEMALE.SAIGE')
            if male_test not in res_files and female_test not in res_files:
                res_files.append(f)
            
    return [os.path.join(res_path,f) for f in res_files]


def return_phenos(pheno_paths):

    phenos = []
    with open(pheno_paths,'rt') as i:
        for line in i:
            pheno = line.strip().split('/')[-1].split('.')[0]
            phenos.append(pheno)

    return phenos

def get_hits(res_path,pheno_paths,g2v_path,out_path,cutoff = 7,lof = 'most_severe',test = True,sep ='\t',blacklist = r2_blacklist):

    make_sure_path_exists(out_path)
    # list of files
    file_list = merge_gender(res_path)
    phenos = return_phenos(pheno_paths)
    # gene to variant dict
    with open(g2v_path,'rb') as i: g2v = pickle.load(i)
    out_file = os.path.join(out_path,lof + '_best_hits.txt')
    tmp_file = os.path.join(out_path,lof + '_best_hits.tmp')
    pretty_print('merging')
    final_phenos =os.path.join(out_path,lof + '_final_phenos.txt')
    rej_phenos = os.path.join(out_path,lof + '_rejected_phenos.txt')
    final = open(final_phenos,'wt')
    rej = open(rej_phenos,'wt')
    with open(tmp_file,'wt') as o:
        if test :
            file_list = file_list[:10]
        #loop through files
        for i,f in enumerate(file_list):
            progressBar(i,len(file_list))
            #get phenotype
            pheno = f.split('.')[0].split('-')[1]
            # check if pheno is blacklisted
            blacklist_check = True
            for bl in blacklist:
                if pheno.startswith(bl):
                    blacklist_check = False
                    rej.write(pheno + ' blacklisted' + '\n')

            #check if pheno in pheno_list
            pheno_check = True if pheno in phenos else False
            if not pheno_check:
                rej.write(pheno + ' missing' + '\n')
                
                
            if blacklist_check and pheno_check:
                final.write(pheno + '\n')
                #file iterator
                iterator = basic_iterator(f,separator = ' ')
                #get releavant info
                header = next(iterator)
                pval_index = header.index('p.value')
                gene_index = header.index('GENE')
                cases_index = header.index('N.Cases')
                controls_index = header.index('N.Controls')
                af_cases_index = header.index('AF.Cases')
                af_controls_index = header.index('AF.Controls')
              
                header[cases_index] = 'ref.count.cases'
                header[controls_index] = 'alt.count.cases'
                header[af_cases_index] = 'ref.count.ctrls'
                header[af_controls_index] = 'alt.count.ctrls'
                
                #loop through data
                for entry in iterator:
                    alt_count_cases = int(round(2*int(entry[cases_index])*float(entry[af_cases_index]),2))
                    ref_count_cases = int(entry[cases_index]) - int(alt_count_cases)                                       
                    assert (alt_count_cases + ref_count_cases == int(entry[cases_index]))
                    
                    alt_count_controls = int(round(2*int(entry[controls_index])*float(entry[af_controls_index]),2))
                    ref_count_controls = int(entry[controls_index]) - int(alt_count_controls)                                       
                    assert (alt_count_controls + ref_count_controls == int(entry[controls_index]))
                    
                    

                    entry[cases_index] = str(ref_count_cases)               
                    entry[controls_index] = str(alt_count_cases)

                    entry[af_cases_index] = str(ref_count_controls)                               
                    entry[af_controls_index] = str(alt_count_controls)

                    gene = entry[gene_index]
                    variants = g2v[gene]
                   
                    entry.append(','.join(variants))
                    out_entry = sep.join([pheno] + entry) + '\n'
                   # print(out_entry)
                    o.write(out_entry)

                
    rej.close()
    final.close()

    final_header = ['pheno'] + header + ['variants']
    pretty_print('sorting')
    with open(out_file,'wt') as o:o.write(sep.join(final_header) + '\n')
        
    cmd = 'sort -g -k ' + str(final_header.index('p.value')+1)  + pad(tmp_file) + ' >> ' + pad(out_file)
    print(cmd)
    tmp_bash(cmd)
    pretty_print('zipping')
    cmd = 'gzip -f ' + out_file
    call(shlex.split(cmd))
    
    
    
def best_hits(resPath ,iPath,lofString = 'hc_lof'):


    lofPath = resPath + '/' + lofString+ '/' 
    header = lofPath + 'columns.txt'
    qqPath = lofPath  +lofString + '_qq_data.txt'
    oPath =  lofPath  +lofString + '_hits.txt'
    resPath += '/' + lofString +'/results/'
    print('fetching data from ' + resPath)
    files = get_filepaths(resPath)
    lines = []
    for f in files:
        pheno = f.split('/')[-1].split('-')[0]
        print(pheno)
        pheno = [pheno]
        with gzip.open(f,'rt') as i:
            for line in i:
                line = line.strip().split('\t')
                pval = line[1]
                gene = line[0]
                variants = g2v[gene]
                try:
                    p = np.float128(pval)
                    pExp = -np.log10(p)
                    line[1] = p
                    lines.append(pheno + line + variants)

                except:
                    pass
    print('sorting...')
    lines = sorted(lines,key = lambda x:x[2])
    with open(header,'rt') as i:
        head = i.readlines()[0]

    head = 'pheno gene ' + head.strip() + ' variants'
    print(head)
    with open(oPath,'wt') as o:
        # add header
        head = head.replace(' ','\t')
        o.write(head+'\n')
        for s in lines:
            oString = '\t'.join([str(elem) for elem in s])
            o.write(oString + '\n')
    cmd = 'gzip -f ' + oPath
    call(shlex.split(cmd))


def get_filepaths(mypath):
    from os import listdir
    from os.path import isfile, join
    return [f for f in listdir(mypath) if isfile(join(mypath, f))]
