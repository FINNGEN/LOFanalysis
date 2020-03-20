#!/usr/bin/env python3

from file_utils import tmp_bash,make_sure_path_exists,file_exists,basic_iterator,return_header
import argparse,json,os
from tempfile import NamedTemporaryFile
   
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
    
    # MERGE RESULTS WITHOUT HEADER
    print('merging results...')
    
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

    print('sorting by pval..')
    sort_col = 1 + header.index('p.value')
    cmd = f' sort -gk {sort_col} {tmp_file} >> {out_file}'
    print(cmd)
    tmp_bash(cmd)
   
    print('zipping...')
    cmd = f"bgzip -f {out_file}"
    print(cmd)
    tmp_bash(cmd)
    
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
        
