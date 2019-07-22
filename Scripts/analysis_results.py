#!/usr/bin/env python3
from Tools.utils import tmp_bash,make_sure_path_exists,get_filepaths,progressBar
import argparse
import os

def download_dict(args):

    args.dict_path = os.path.join(args.outpath,'dict')
    make_sure_path_exists(args.dict_path)
    cmd = f"gsutil cp gs://fg-cromwell/{args.workflow}/{args.id}/call-gene_matrix/**/*dict.txt {args.dict_path}"
    print(cmd)
    tmp_bash(cmd)

    cmd = f"gsutil cp gs://fg-cromwell/{args.workflow}/{args.id}/call-merge_matrix/**/*matrix.tsv {os.path.join(args.outpath,'gene_matrix.tsv')}"
    print(cmd)
    tmp_bash(cmd)
    
    cmd = f"jq -s  'reduce .[] as $item ({{}}; . * $item)' {args.dict_path}/* > {os.path.join(args.outpath,'gene_dict.json')}"
    print(cmd)
    tmp_bash(cmd)

def download_saige(args):
    cmd = f"gsutil -m cp gs://fg-cromwell/{args.workflow}/{args.id}/call-pheno_saige/**/*SAIGE.txt {args.saige_path}"
    tmp_bash(cmd)
    

def saige_merge(args):
    
    files = get_filepaths(args.saige_path)
    out_file = os.path.join(args.outpath,args.lof + "_gene_results.txt")
    tmp_file = os.path.join(args.outpath,args.lof + "_tmp.txt")

    if os.path.isfile(out_file):
        return

    else:
        print('merging results...')
        with open(tmp_file,'wt') as o:
            for i,saige_file in enumerate(files):
                progressBar(i,len(files))
                pheno = os.path.basename(saige_file).split('.')[0]
                with open(saige_file) as i:
                    header = next(i)
                    for line in i:
                        o.write(pheno +" " +line)

        with open(out_file,'wt') as o:
            o.write("PHENO " +header)

        print('sorting...')
        pval_index = header.split(' ').index('p.value') + 2
        cmd = f"sort -gk {pval_index} {tmp_file}  >> {out_file}"
        tmp_bash(cmd)
    
    
if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description="Deal with lof results")

    parser.add_argument('--outpath', type=str, help='output path',required = True)
    parser.add_argument('--id', type=str, help='worflow id',required = True)
    parser.add_argument('--workflow', type=str, help='worflow name',default = "LOF_saige")
    
    subparsers = parser.add_subparsers(help='help for subcommand',dest ="command")
    dict_parser = subparsers.add_parser('dict', help='download gene dictionaries and merge them')
    saige_download_parser = subparsers.add_parser('saige_download', help='download saige results')
    saige_merge_parser = subparsers.add_parser('saige', help='download saige results')
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
