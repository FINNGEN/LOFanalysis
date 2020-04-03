#!/usr/bin/env python3

from file_utils import tmp_bash,make_sure_path_exists,file_exists,basic_iterator,return_header,get_filepaths,mapcount,StoreDictKeyPair
import argparse,json,os,scipy,subprocess,shlex,shutil
from tempfile import NamedTemporaryFile
from pathlib import Path
import numpy as np
from collections import defaultdict,Counter
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(palette='Reds')

# fields that need to be there
required_fields =['PHENO', 'GENE', 'p.value', 'BETA', 'SE', 'AC', 'AF', 'N', 'N.Cases', 'N.Controls', 'AF.Cases', 'AF.Controls']
# fields that the script adds
additional_fields = ['ref.count.cases','alt.count.cases','ref.count.controls','alt.count.controls','variants']

def saige_merge(args):

   
    tmp_file = os.path.join(args.outpath, "tmp.txt")
    out_file = os.path.join(args.outpath, f"{args.prefix}_gene_results.txt")       
    args.out_file = out_file+'.gz'
    
    args.json = os.path.join(args.outpath, f"{args.prefix}_gene.json")       
    subprocess.call(shlex.split(f"cp {args.dict_file} {args.json}"))

    print(args.out_file)
    if os.path.isfile(args.out_file):
        print('all files already generated')
        return

    sanitize_inputs(args)
    print("editing metadata...",args.out_file)

    # MERGE JSON FILES
    with open(args.json) as infile: g_dict = json.load(infile)
    j = {}
    for gene in g_dict:
        j[gene] =  ','.join(g_dict[gene])
    
    saige_iterator = basic_iterator(args.saige_file)
    with open(tmp_file,'w') as t,open(out_file,'w') as o :
        # gather all required column indexes, replacing values if needed
        header = next(saige_iterator)
        for key in args.fields_dict:
            header[header.index(key)] = args.fields_dict[key]          
        # write new header adding new columns
        new_header = header+ additional_fields
        o.write('\t'.join(new_header) + '\n')
        # relevant fields that i need to process while parsing the file
        gene_idx = new_header.index("GENE")
        idx_list = [new_header.index(elem) for elem in ['N.Cases','N.Controls','AF.Cases','AF.Controls']]
        for line in saige_iterator:
            gene = line[gene_idx]
            cases,controls,af_cases,af_controls = list(map(float,[line[idx] for idx in idx_list]))
            # create new required entries (additional fields)
            alt_count_cases = af_cases * cases*2
            ref_count_cases = cases - alt_count_cases
            alt_count_controls = af_controls*controls*2
            ref_count_controls = controls - alt_count_controls
            res = list(map(int,[ref_count_cases,alt_count_cases,ref_count_controls,alt_count_controls]))
            # create out_line
            new_line = line +  list(map(str,res)) + [j[gene]] 
            assert len(new_line) == len(new_header)
            t.write('\t'.join(new_line) + '\n')

            
    sort_col = 1 + header.index('p.value')
    cmd = f'sort -gk {sort_col} {tmp_file} >> {out_file}'
    print('sorting by pval:',cmd)
    tmp_bash(cmd)
   

    cmd = f"bgzip -f {out_file}"
    print('zipping:',cmd)
    tmp_bash(cmd)



##############################
#----------PLOTTING----------#
##############################

def return_pheno_counts(f,bin_col):

    #create dict that maps pheno to number of cases
    header = return_header(f)
    count_iterator = basic_iterator(f,columns = [header.index(elem) for elem in ['PHENO','p.value',bin_col]],skiprows = 1)
    cases_count,pheno_pvals = {},defaultdict(list)
    for pheno,pval,cases in count_iterator:
        cases_count[pheno] = int(float(cases))
        pheno_pvals[pheno].append(float(pval))
    return cases_count,pheno_pvals  

def plot_qq(f,cases_count = None,pheno_pvals = None,bin_col = 'N.Cases'):

    args.qq = f.replace('results.txt.gz','qq.pdf')
    print(args.qq)
    if os.path.isfile(args.qq):
        print('qq plot already generated')
        return
    
    # count cases for each pheno
    if not cases_count and not pheno_pvals:
        cases_count,pheno_pvals = return_pheno_counts(f,bin_col)
        
    with open(os.path.join(args.outpath,'phenos.txt'),'wt') as o:
        for pheno in pheno_pvals:o.write(pheno + '\n')
        
    bins = [1,200,500,1000,10000,int(max(cases_count.values()))+1]
    labels = [str((bins[i],bins[i+1])) for i in range(len(bins)-1)]

    # create bins 
    print('binning phenotypes')
    pheno_bin = {} 
    for pheno in cases_count:
        pheno_bin[pheno] = labels[np.digitize(cases_count[pheno],bins)-1] 
    
    print('binning pvalues')
    pval_bin =defaultdict(list)
    for pheno in pheno_pvals:
        pval_bin[pheno_bin[pheno]] += pheno_pvals[pheno]
        
    label_count = Counter(pheno_bin.values())
    print(label_count)
    
    print('plotting...')
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1,1) 
    ax = fig.add_subplot(gs[0,0])

    max_list = []
    for label in labels:
        print(label)
        pvalues = pval_bin[label]
        n = len(pvalues)
        if n:
            pvalues = np.array(pvalues)
            inflation_factor = lambda_inflation(pvalues)
            pvalues = sorted(-np.log10(pvalues))
            theory_p = sorted(-np.log10((np.arange(n)+0.5) / n))
            max_list.append(max(theory_p))
            ax.plot(theory_p,pvalues,label = label + f" | {label_count[label]} | {inflation_factor}" )

    # plot diagonal
    diag = np.linspace(0,max_list,100)
    ax.plot(diag,diag,'--',color = 'k',linewidth=0.5)

    plt.title("QQ plot for gene based saige")
    ax.set_xlabel(r'expected -$log_{10}(p)$')
    ax.set_ylabel(r'observed -$log_{10}(p)$') 
    ax.set_xlim((0,max(max_list)))
    ax.set_ylim((0,max(max_list)))
    
    leg = plt.legend(loc = 'lower right',fontsize = 6)
    leg.set_title('cases range | pheno count | lambda',prop={'size':7})
    
    fig.savefig(args.qq)
    plt.close()

def lambda_inflation(pvalues):
    stats = scipy.stats.chi2.ppf(1-pvalues,1)
    inflation_factor = round(np.median(stats) /scipy.stats.chi2.ppf(0.5, 1) ,2)
    return inflation_factor

def release(args):

    doc_path = os.path.join(args.outpath,'documentation')
    data_path = os.path.join(args.outpath,'data')
    for path in [doc_path,data_path]:
        make_sure_path_exists(path)
        for f in get_filepaths(path): os.remove(f) # clean path else shutil.copy might fail
    for f in [args.out_file,args.json]:
        shutil.copy(f,os.path.join(data_path,os.path.basename(f)))
    for f in [args.qq]:
        shutil.copy(f,os.path.join(doc_path,os.path.basename(f)))
        

    # get variant/gene counts & summary
    with open(args.json) as infile: g_dict = json.load(infile)
    genes = len(g_dict.keys())
    variants = np.sum([len(g_dict[k]) for k in g_dict])
    
    gene_summary = "|# of variants in gene | count|\n|--|--|\n"
    counter = Counter([len(g_dict[k]) for k in g_dict])
    for n in sorted(counter.keys(),reverse=True):
        gene_summary += "|" + "|".join(list(map(str,[n,counter[n]]))) + '|\n'

    phenos = mapcount(os.path.join(args.outpath,'phenos.txt'))
    
    parent_path = Path(os.path.realpath(__file__)).parent.parent
    readme = os.path.join(parent_path,'Data','lof.README')    
    with open( os.path.join(args.outpath,args.prefix + '_lof_readme'),'wt') as o, open(readme,'rt') as i:
        word_map = {
            '[PREFIX]':args.prefix,
            '[N_VARIANTS]' : variants,
            '[N_GENES]' : genes,
            '[GENE_SUMMARY]':gene_summary,
            '[N_PHENOS]': phenos
            }
        
        word_map.update(args.release_dict)
       
        for line in i:
            for kw in word_map:
                if kw in line:
                    line = line.replace(kw,str(word_map[kw]))
                    
            o.write(line)
            
def sanitize_inputs(args):
    # the values are the required entries on pheweb
    args.fields_dict = {}
    for i,elem in enumerate([args.pheno_col,args.gene_col,args.pval_col,args.beta_col,args.se_col,args.ac_col,args.af_col,args.n_col,args.cases_col,args.controls_col,args.af_cases_col,args.af_controls_col]):
        args.fields_dict[elem] = required_fields[i]
    
    print(args.fields_dict.values())

    header = return_header(args.saige_file)

    missing_keys = [key for key in args.fields_dict if key not in header]
    if missing_keys:
        print(header)
        raise ValueError(f"{missing_keys} missing in header")
      
    print('all columns in header')

if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description="Deal with lof results")

    parser.add_argument('-o','--outpath', type=str, help='output path',default = os.getcwd())      
    parser.add_argument('--saige-file', type = file_exists)
    parser.add_argument('--dict-file', type = file_exists)
    parser.add_argument('--prefix', type = str, default = "finngen")
    parser.add_argument('--release', action = 'store_true')
    #COLUMN NAMES
    parser.add_argument('--pheno-col', type = str, default = "PHENO")
    parser.add_argument('--gene-col', type = str, default = "GENE")
    parser.add_argument('--ac-col', type = str, default = "AC_Allele2")
    parser.add_argument('--af-col', type = str, default = "AF_Allele2")
    parser.add_argument('--n-col', type = str, default = "N")
    parser.add_argument('--beta-col', type = str, default = "BETA")
    parser.add_argument('--se-col', type = str, default = "SE")
    parser.add_argument('--pval-col', type = str, default = "p.value")    
    parser.add_argument('--cases-col', type = str, default = "N.Cases")
    parser.add_argument('--controls-col', type = str, default = "N.Controls")
    parser.add_argument('--af-cases-col', type = str, default = "AF.Cases")
    parser.add_argument('--af-controls-col', type = str, default = "AF.Controls")
    parser.add_argument("--release_dict", action=StoreDictKeyPair, metavar="KEY1=VAL1,KEY2=VAL2...")

    args = parser.parse_args()
    
 
    saige_merge(args)
    plot_qq(args.out_file)
    if args.release:
        print('release!')
        release(args)
