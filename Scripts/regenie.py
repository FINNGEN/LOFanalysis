from utils_regenie import get_exitcode_stdout_stderr,get_param_list,check_error_regenie,read_phenos
from itertools import product
from file_utils import progressBar,file_exists,make_sure_path_exists,mapcount,get_filepaths
import multiprocessing,time,argparse,os
cpus = multiprocessing.cpu_count()

# DEFAULTS FOR TIME BEING
covar_list=" SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC{1:10},IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm,BATCH_DS18_MIGRAINE_1_norm,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm,BATCH_DS21_SUPER_2_norm,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm "


def regenie_cmd(pred,pgen,set_list,mask,annotation,covar_file, aaf_bins,aaf_file, covars = covar_list, regenie_args = " --bsize 200   --aaf-bins 0.01,0.1,0.5"):
    """out
    Produces regenie run commands
    """
    basic_cmd = f"regenie --step 2 --bt --ref-first --firth --approx --threads 1  "

    # regenie burden required files
    aaf = " --aaf-file " + aaf_file if aaf_file else ""
    bins = " --aaf-bins " + " ".join(map(str,aaf_bins)) if aaf_bins else ""
    basic_cmd += f" --pred {pred} --pgen {pgen.replace('.pgen','')} --mask-def {mask} --set-list {set_list} --anno-file {annotation} "
    basic_cmd += f" --covarFile {covar_file} --covarColList {covars} --phenoFile {covar_file} "
    basic_cmd += f" {aaf} {bins} {regenie_args} "
    print(basic_cmd)
    return basic_cmd


def regenie_run(basic_cmd,args,pheno,chrom):
    """
    Single regenie run, it will retry until it successfully runs, else it will log
    """

    # update command with chrom/pheno info
    out_file = os.path.join(args.out,args.name)
    cmd = basic_cmd + f" --chr {chrom} --phenoColList {pheno} --out {out_file}_{chrom}"
    # re-run fun until either success or exceeded max retries
    retries,success = 0,False
    while retries < args.max_retries and not success:
        exitcode,stdout,err = get_exitcode_stdout_stderr(cmd)
        if check_error_regenie(stdout):retries += 1
        else:success = True

    print(f"\n{pheno} {chrom} Success:{success} ")

    # log results
    with open(f"{out_file}_{chrom}_{pheno}.log",'wt') as o,open(out_file + "_failed.txt",'at') as f:
        o.write(stdout)
        if not success:
            f.write(f"{chrom}_{pheno}\n")

    return pheno,chrom,success,stdout

def wrapper_func(args):
    return regenie_run(*args)



def merge_results(args):
    params = product(args.phenos,args.chrom)
    out_regenie = os.path.join(args.out,f"{args.name}.regenie")
    out_header = os.path.join(args.out,f"{args.name}.header")
    out_log = os.path.join(args.out,f"{args.name}.log")
    print(out_regenie)
    with open(out_regenie,'wt') as o,open(out_header,'wt') as h,open(out_log,'wt') as log:
        for pheno,chrom in params:
            result = os.path.join(args.out,f"{args.name}_{chrom}_{pheno}.regenie")
            log_file =  os.path.join(args.out,f"{args.name}_{chrom}_{pheno}.log")
            with open(result) as i:
                # dump header
                header = [next(i) for j in range(2)]
                header[1] = "PHENO " + header[1]
                for line in i:o.write(f"{pheno} {line}")
            with open(log_file )as i:
                for line in i:
                    log.write(line)
        for elem in header:h.write(elem)

    return

def main(args):

    # clear previous failed runs
    os.remove(args.out + f"{args.name}_failed.txt") if os.path.isfile(args.out + f"{args.name}_failed.txt") else None

    basic_cmd =regenie_cmd(args.pred ,args.pgen,args.sets,args.mask , args.annotation,args.pheno_file,args.aaf_bins,args.aaf_file,args.covariates,args.regenie_args)
    params = [(basic_cmd,args,pheno,chrom) for chrom,pheno in get_param_list(args.out,args.name,args.pred,args.chrom,args.phenos,args.force)]

    if not params:
        print("All phenos ran")
        return

    print(f"{len(params)} combinations to loop over.")
    pools = multiprocessing.Pool(cpus - 1)
    results = pools.map_async(wrapper_func,params)
    while not results.ready():
        time.sleep(10)
        progressBar(len(params) - results._number_left,len(params))

    progressBar(len(params) - results._number_left,len(params))
    pools.close()





if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Regenie LOF burden")
    # required args
    parser.add_argument('-o',"--out",type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--pgen",type = file_exists,help = "Pgen file",required = True)
    parser.add_argument("--pheno-file",type = file_exists,help = "Regenie pheno file",required = True)
    parser.add_argument("--mask",type = file_exists,help = "Regenie mask file",required = True)
    parser.add_argument("--annotation",type = file_exists,help = "Regenie annotation file",required = True)
    parser.add_argument("--sets",type = file_exists,help = "Regenie sets file",required = True)
    parser.add_argument("--pred",type = file_exists,help = "Regenie pred file",required = True)
    parser.add_argument("--aaf-file",type = file_exists,help = "Regenie aaf file",required = False)
    parser.add_argument("--aaf-bins",nargs="+",type=float,required = False)
    parser.add_argument("--covariates",type = str, help = "Covariates",default = covar_list)
    parser.add_argument('-c','--chrom',nargs="+",type=str,choices = list(map(str,range(1,24))) + ['X'],help = 'chromosome number',  default = range(1,23))
    parser.add_argument("--regenie-args",type = str,help = "Extra kwarg to pass to regenie",default = " --bsize 200   " )
    parser.add_argument("--name",type =str, help = "output root",default = "finngen_lof")
    parser.add_argument('--cpus',type = int,help = 'number of parallel processes to run', default = cpus)
    parser.add_argument('--max-retries',type = int,help = 'Maximum regenie retries', default = 2)
    parser.add_argument('--test',action = 'store_true',help = 'Runs test version')
    parser.add_argument('--force',action = 'store_true',help = 'Re runs everything')

    args=parser.parse_args()
    make_sure_path_exists(args.out)

    #print(args)
    if args.test:
        args.chrom = ["21","22"]

    args.phenos = read_phenos(args.pred)

    print(args)
    main(args)
    merge_results(args)
