#!/usr/bin/env python3

from file_utils import tmp_bash,make_sure_path_exists,file_exists,basic_iterator,return_header
import argparse,json,os,scipy
from tempfile import NamedTemporaryFile
import numpy as np
from collections import defaultdict,Counter
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
#sns.set(palette='Set2')
sns.set(palette='Reds')

# fields that need to be there
required_fields =['PHENO', 'GENE', 'p.value', 'BETA', 'SE', 'AC', 'AF', 'N', 'N.Cases', 'N.Controls', 'AF.Cases', 'AF.Controls']
# fields that the script adds
additional_fields = ['ref.count.cases','alt.count.cases','ref.count.controls','alt.count.controls','variants']

def saige_merge(args):
   
    # MERGE JSON FILES
    with open(args.dict_file) as infile: g_dict = json.load(infile)
    j = {}
    for gene in g_dict:
        j[gene] =  ','.join(g_dict[gene])
    

    print("editing metadata...")
    tmp_file = os.path.join(args.outpath, "tmp.txt")
    out_file = os.path.join(args.outpath, f"{args.prefix}_gene_results.txt")       
    
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

    out_file += '.gz'
    plot_qq(out_file,savefig = out_file.replace('.txt.gz','.pdf'))



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
    


def plot_qq(f,cases_count = None,pheno_pvals = None,savefig = '/mnt/disks/r5/lof/results/test.pdf',bin_col = 'N.Cases'):

    # count cases for each pheno
    if not cases_count and not pheno_pvals:
        cases_count,pheno_pvals = return_pheno_counts(f,bin_col)

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
    
    fig.savefig(savefig)
    plt.close()
    print(savefig,'done')      

def lambda_inflation(pvalues):
    stats = scipy.stats.chi2.ppf(1-pvalues,1)
    inflation_factor = round(np.median(stats) /scipy.stats.chi2.ppf(0.5, 1) ,2)
    return inflation_factor




if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description="Deal with lof results")

    parser.add_argument('-o','--outpath', type=str, help='output path',default = os.getcwd())      
    parser.add_argument('--saige-file', type = file_exists)
    parser.add_argument('--dict-file', type = file_exists)
    parser.add_argument('--prefix', type = str, default = "finngen")
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

    args = parser.parse_args()

    # the values are the required entries on pheweb
    args.fields_dict = {}
    for i,elem in enumerate([args.pheno_col,args.gene_col,args.pval_col,args.beta_col,args.se_col,args.ac_col,args.af_col,args.n_col,args.cases_col,args.controls_col,args.af_cases_col,args.af_controls_col]):
        args.fields_dict[elem] = required_fields[i]
    
    print(args.fields_dict.values())

    header = return_header(args.saige_file)
    
    for key in args.fields_dict:
        if key not in header:
            raise ValueError(f"{key} not found in file header")
      
    print('all columns in header')
    saige_merge(args)    
        
