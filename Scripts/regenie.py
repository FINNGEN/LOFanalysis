from utils_regenie import get_exitcode_stdout_stderr
from file_utils import progressBar,file_exists,make_sure_path_exists
import multiprocessing,time,argparse,os
from itertools import product
from functools import partial
cpus = multiprocessing.cpu_count()

# DEFAULTS FOR TIME BEING
covar_list=" SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC{1:10},IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm,BATCH_DS18_MIGRAINE_1_norm,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm,BATCH_DS21_SUPER_2_norm,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm "


def regenie_cmd(pred,pgen,set_list,mask,annotation,covar_file, out, covars = covar_list, regenie_args = " --bsize 200   --aaf-bins 0.01,0.1,0.5" ,chrom_list = range(1,23)):
    """
    Produces regenie run commands
    """
    basic_cmd = f"regenie --step 2 --bt --ref-first --firth --approx --threads 1  --out {out}_CHROM"

    # regenie burden required files
    basic_cmd += f" --pred {pred} --pgen {pgen.replace('.pgen','')} --mask-def {mask} --set-list {set_list} --anno-file {annotation}"
    basic_cmd += f" --covarFile {covar_file} --covarColList {covars} --phenoFile {covar_file} "
    basic_cmd += f" {regenie_args} "
    print(basic_cmd)

    return basic_cmd


def regenie_run(basic_cmd,pheno,chr,max_retries = 10):
    """
    Single regenie run, it will retry until it successfully runs, else it will log
    """

    # update command with chrom/pheno info
    cmd = basic_cmd + f" --chr {chr} --phenoColList {pheno} "
    cmd = cmd.replace("CHROM",str(chr))

    # re-run fun until either success or exceeded max retries
    retries,success = 0,False
    while retries < max_retries and not success:
        exitcode,stdout,err = get_exitcode_stdout_stderr(cmd)

        if check_error_regenie(stdout):
            retries += 1
        else:
            success = True

    print(f"\nSuccess:{success} {pheno} {chr}")

    return pheno,chr,success,stdout

def wrapper_func(args):
    return regenie_run(*args)

def check_error_regenie(stdout):
    return "ERROR" in stdout

def read_phenos(pred):
    with open(pred) as i: phenos = [elem.strip().split()[0] for elem in i]
    return phenos


def log_results(out,res):
    """
    Write stdout of each run to log and stores failed runs
    """
    f = open(out + "_failed.txt",'wt')
    for entry in res.get():
        pheno,chrom,success,stdout = entry
        with open(f"{out}_{pheno}_{chrom}.log",'wt') as o:
            o.write(stdout)
        if not success:
            f.write(f"{pheno}_{chr}\n")
    f.close()

def main(pred,pgen,set_list,mask,annotation,covar_file, out, covar_list = covar_list,regenie_args = " --bsize 200   --aaf-bins 0.01,0.1,0.5" ,chrom_list = range(1,23),max_retries = 10):

    basic_cmd =regenie_cmd(pred ,pgen,set_list,mask , annotation,covar_file,out,covar_list,regenie_args,chrom_list)
    phenos = read_phenos(pred)

    params = list(product([basic_cmd],phenos,chrom_list,[max_retries]))

    pools = multiprocessing.Pool(cpus - 1)
    results = pools.map_async(wrapper_func,params)
    while not results.ready():
        time.sleep(10)
        progressBar(len(params) - results._number_left,len(params))

    progressBar(len(params) - results._number_left,len(params))
    pools.close()
    log_results(out,results)


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

    parser.add_argument("--covariates",type = str, help = "Covariates",default = covar_list)
    parser.add_argument('-c','--chrom',choices = list(map(str,range(1,24))) + ['X'],help = 'chromosome number',  default = range(1,22))
    parser.add_argument("--regenie-args",type = str,help = "Extra kwarg to pass to regenie",default = " --bsize 200   --aaf-bins 0.01,0.1,0.5" )
    parser.add_argument("--name",type =str, help = "output root",default = "finngen_lof")
    parser.add_argument('--cpus',type = int,help = 'number of parallel processes to run', default = cpus)
    parser.add_argument('--max-retries',type = int,help = 'Maximum regenie retries', default = 10)
    parser.add_argument('--test',action = 'store_true',help = 'Runs test version')

    args=parser.parse_args()
    make_sure_path_exists(args.out)
    args.out = os.path.join(args.out,args.name)
    #print(args)
    if args.test:
        args.chrom = [22]

    print(args)
    main(args.pred,args.pgen,args.sets,args.mask,args.annotation,args.pheno_file,args.out,args.covariates,args.regenie_args,args.chrom,args.max_retries)
