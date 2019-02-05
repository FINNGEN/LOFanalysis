from file_utils import *

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


def get_hits(res_path,g2v_path,out_path,cutoff = 7,lof = 'most_severe',test = True,sep =','):

    make_sure_path_exists(out_path)
    # list of files
    file_list = merge_gender(res_path)
    # gene to variant dict
    with open(g2v_path,'rb') as i: g2v = pickle.load(i)
    out_file = os.path.join(out_path,lof + '_best_hits.txt')
    tmp_file = os.path.join(out_path,lof + '_best_hits.tmp')
    pretty_print('merging')
    with open(tmp_file,'wt') as o:
        if test :
            file_list = file_list[:10]
        #loop through files
        for i,f in enumerate(file_list):
            progressBar(i,len(file_list))
            #get phenotype
            pheno = f.split('.')[0].split('-')[1]
            #file iterator
            iterator = basic_iterator(f,separator = ' ')
            #get releavant info
            header = next(iterator)
            pval_index = header.index('p.value')
            gene_index = header.index('GENE')
         
            
            #loop through data
            for entry in iterator:              
                gene = entry[gene_index]
                variants = g2v[gene]
                entry.append(' '.join(variants))
                out_entry = sep.join([pheno] + entry) + '\n'
                if cutoff:
                    pval = np.float128(entry[pval_index])
                    exp = -np.log10(pval)
                    if exp > cutoff:
                        o.write(out_entry)
                else:
                    o.write(out_entry)
                    
    
    final_header = ['pheno'] + header + ['variants']
    pretty_print('sorting')
    with open(out_file,'wt') as o:o.write(sep.join(final_header) + '\n')
    cmd = 'sort -g -k ' + str(final_header.index('p.value')+1)  + ' -t' + pad(sep) + pad(tmp_file) + ' >> ' + pad(out_file)
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
