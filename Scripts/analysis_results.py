#!/usr/bin/env python3

from Tools.utils import tmp_bash,make_sure_path_exists,get_filepaths,progressBar,file_exists,basic_iterator
import argparse
import json
import os
from tempfile import NamedTemporaryFile

def download_dict(args):

    args.dict_path = os.path.join(args.outpath,'dict')
    make_sure_path_exists(args.dict_path)
    cmd = f"gsutil cp gs://fg-cromwell/{args.workflow}/{args.id}/call-gene_matrix/**/*dict.txt {args.dict_path}"
    print(cmd)
    tmp_bash(cmd)

    cmd = f"gsutil cp gs://fg-cromwell/{args.workflow}/{args.id}/call-merge_matrix/**/*matrix.tsv {os.path.join(args.outpath,'gene_matrix.tsv')}"
    print(cmd)
    tmp_bash(cmd)
    
    #cmd = f"jq -s  'reduce .[] as $item ({{}}; . * $item)' {args.dict_path}/* > {os.path.join(args.outpath,'gene_dict.json')}"
    #print(cmd)
    #tmp_bash(cmd)

def download_saige(args):
    cmd = f"gsutil -m cp gs://fg-cromwell/{args.workflow}/{args.id}/call-pheno_saige/**/*SAIGE.txt {args.saige_path}"
    tmp_bash(cmd)
    

def saige_merge(args):

    # LOAD DICT & SAIGE RESULTS
    if not args.saige_files:
        args.saige_files = [f for f in get_filepaths(args.saige_path) if f.endswith('SAIGE.txt')]
    else:
        with open(args.saige_files) as i :
            args.saige_files = [elem.strip() for elem in i.readlines()]
    print(f"{len(args.saige_files)} saige files.")
    if not args.dict_files:
        dict_path = os.path.join(args.outpath,'dict')
        args.dict_files = [f for f in get_filepaths(dict_path) if f.endswith('dict.txt')]
    else:
        with open(args.dict_files) as i :
            args.dict_files = [elem.strip() for elem in i.readlines()]
    print(f"{len(args.dict_files)} dict files.")
    
    # MERGE JSON FILES
    g = {}
    for f in args.dict_files:
        print(f)
        with open(f) as infile:
            g.update(json.load(infile))

    gene_json = os.path.join(args.outpath,args.lof + "_gene_dict.json")
    with open(gene_json, "w") as outfile:
        json.dump(g, outfile)
    with open(gene_json) as i :g_dict = json.load(i)

    # MERGE RESULTS WITHOUT HEADER
    print('merging results...')                
    tmp_fix = NamedTemporaryFile(delete=True)
    with open(tmp_fix.name,'w') as o:
        for i,saige_file in enumerate(args.saige_files):
            progressBar(i,len(args.saige_files))
            pheno = os.path.basename(saige_file).split('.')[0]
            with open(saige_file) as i:
                header = next(i)
                for line in i:
                    o.write(pheno +" " +line)

    # CREATE TMP FILE & SORT BY PVAL
    tmp_file = os.path.join(args.outpath,args.lof + "_tmp.txt")
    with open(tmp_file,'wt') as t:
        t.write("\t".join(["PHENO"] + header.strip().split(" ")) + '\n')
    print('sorting...')
    pval_index = header.split(' ').index('p.value') + 2
    cmd = f"sort -gk {pval_index} {tmp_fix.name} | tr [:blank:] \\\\t >> {tmp_file}"
    print(cmd)
    tmp_bash(cmd)

    # ADD METADATA TO FINAL FILE
    out_file = os.path.join(args.outpath,args.lof + "_gene_results.txt")   
    with open(out_file,'wt') as o :
        res_iterator = basic_iterator(tmp_file)
        header = next(res_iterator)
        # replace column headers
        idx_list =[ header.index(elem) for elem in ['N.Cases','N.Controls','AF.Cases','AF.Controls']]
        for i,elem in enumerate(['ref.count.cases','alt.count.cases','ref.count.controls','alt.count.controls']):
            header[idx_list[i]] = elem
          
        header.append('variants')
        o.write('\t'.join(header) +'\n')

        count = 0
        for entry in res_iterator:
            cases,controls,af_cases,af_controls = list(map(float,[entry[ix] for ix in idx_list]))
            alt_cases = af_cases * cases*2
            ref_cases = cases - alt_cases
            alt_controls = af_controls*controls*2
            ref_controls = controls - alt_controls
            res = list(map(str,[ref_cases,alt_cases,ref_controls,alt_controls]))
            for i,ix in enumerate(idx_list):
                entry[ix] = res[i]

            variants = ','.join(g_dict[entry[1]])
            entry.append(variants)
            o.write('\t'.join(entry) +'\n')
            
    os.remove(tmp_file)
    print('zipping...')
    cmd = f"bgzip {out_file}"
    tmp_bash(cmd)
    
if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description="Deal with lof results")

    parser.add_argument('-o','--outpath', type=str, help='output path',required = True)
    parser.add_argument('--id', type=str, help='worflow id')
    parser.add_argument('--workflow', type=str, help='worflow name',default = "LOF_saige")
    
    subparsers = parser.add_subparsers(help='help for subcommand',dest ="command")
    dict_parser = subparsers.add_parser('dict', help='download gene dictionaries and merge them')
    saige_download_parser = subparsers.add_parser('saige_download', help='download saige results')

    saige_merge_parser = subparsers.add_parser('saige', help='download saige results')
    saige_merge_parser.add_argument('--saige_files', type = file_exists)
    saige_merge_parser.add_argument('--dict_files', type = file_exists)
    saige_merge_parser.add_argument('--lof', type = str, default = "most_severe")
    args = parser.parse_args()

    
    if args.command == 'dict':
        download_dict(args)

    if args.command == 'saige_download':
        args.saige_path = os.path.join(args.outpath,'saige')
        make_sure_path_exists(args.saige_path)
        download_saige(args)

    if args.command == 'saige':
        args.saige_path = os.path.join(args.outpath,'saige')
        make_sure_path_exists(args.saige_path)
     

        saige_merge(args)    
