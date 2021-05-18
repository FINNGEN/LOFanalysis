from collections import defaultdict
from subprocess import Popen,PIPE
import shlex,os,gzip
import numpy as np
from file_utils import mapcount
from itertools import product

def copy_nulls_with_name(cromwell_id,dest_bucket = "gs://r7_data/regenie/nulls/",cromwell_machine = "fg-cromwell_fresh",workflow = 'regenie'):


    # here i build the command to retrieve a list of lists of files.
    pred_lists_cmd = f"gsutil ls gs://{cromwell_machine}/{workflow}/{cromwell_id}/call-step1/**/*pred.list"
    print(pred_lists_cmd)
    exitcode, out, err = get_exitcode_stdout_stderr(pred_lists_cmd)
    # now i loop through each list and build a mapping
    for pred_list in [elem for elem in out.split('\n') if elem]:
        print(f"List of nulls {pred_list}")

        shard = pred_list.split("shard-")[1].split('/')[0]
        print(f"shard number :{shard}")

        # now, for this shard, we can map the file id to the pheno
        _,out,_ = get_exitcode_stdout_stderr(f"gsutil cat {pred_list}")
        pheno_mapping = build_pheno_mapping(out)
        print("id-pheno mapping built")

        # now i need to get all the loco.gz files

        # this is the generic command which takes into account globbing and attempts.
        null_lists_cmd = f"gsutil ls gs://{cromwell_machine}/{workflow}/{cromwell_id}/call-step1/shard-{shard}/**/*_ID.loco.gz"
        # for each f_id i return the specific file and move it to dest_bucket
        for f_id in pheno_mapping:
            print(f_id)
            # now i retrieve the loco specific to this id
            _,loco_path,_ = get_exitcode_stdout_stderr(f"{null_lists_cmd.replace('ID',f_id)}")
            # and i proceed to move it
            exitcode,out,err = get_exitcode_stdout_stderr(f"gsutil cp {loco_path} {dest_bucket}{pheno_mapping[f_id]}.loco.gz")
            # check if there's any problem
            if not exitcode: print(err)

    return



def build_pheno_mapping(output):
    """
    Based on the pred_list it maps numerical id to pheno
    """
    pheno_mapping = {}
    for entry in [elem for elem in output.split('\n') if elem]:
        pheno,f = entry.strip().split()
        f_id = f.split('.loco.gz')[0].split('_')[-1]
        pheno_mapping[f_id] = pheno
    return pheno_mapping


def get_exitcode_stdout_stderr(cmd):
    """
    Execute the external command and get its exitcode, stdout and stderr.
    """
    args = shlex.split(cmd)

    proc = Popen(args, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    exitcode = proc.returncode
    #
    return exitcode, out.decode("utf-8").strip() , err



def annotation_filter(annot_file,variant_file,out_root):
    """
    Function that filters the small annotation file to only include variants in the lof pgen.
    """
    variants = set(np.loadtxt(variant_file,dtype = str))
    with gzip.open(annot_file) as i,open(out_root + "annotation.tsv",'wt') as o:
        next(i)
        for line in i:
            variant,_,gene,most_severe = line.strip().decode("utf-8").split()
            variant = "chr" + "_".join(variant.split(':'))
            if variant in variants:
                o.write("\t".join([variant,gene,most_severe]) + "\n")


def filter_annotation(annot_file,tags,out_root = None):
    """
    Returns all variants in genes where at least one variant has consequences in tags
    """
    with open(tags) as f: tag_list = [elem.strip() for elem in f.readlines()]
    print(tag_list)

    if not out_root: out_root = tags.replace('_tags.txt','')


    out_genes = out_root + '_genes.txt'
    print(out_genes)
    if not os.path.isfile(out_genes):
        with open(annot_file) as i:
            genes = []
            # read first line and use data as starting point
            for line in i:
                variant,gene,consequence = line.strip().split()
                if consequence in tag_list:
                    genes.append(gene)

        genes = list(set(genes))
        with open(out_genes,'wt') as o:
            for gene in genes:
                o.write(gene + '\n')
        print('genes saved')
    else:
        with open(out_genes) as i: genes = [elem.strip() for elem in i.readlines()]
        print('genes already saved')


    out_variants = out_root + '_variants.txt'
    with open(annot_file) as i,open(out_variants,'wt') as o:
        # read first line and use data as starting point
        gene_data = gene_iterator(annot_file,genes)
        for data in gene_data:
            assert data[0].strip().split()[1] in genes
            for line in data:
                o.write(line)
    return



def gene_iterator(annot_file,genes):

    with open(annot_file) as i:

        line = next(i)
        variant_data = line.strip().split()
        current_gene = variant_data[1]
        gene_data = [line]

        for line in i:
            # read new data
            variant_data = line.strip().split()
            gene = variant_data[1]
            # if same gene, append data
            if gene == current_gene:
                gene_data.append(line)
            else:
                # we've hit a new gene, time to yield old data and initialize new block
                if current_gene in genes:
                    yield gene_data
                # in any case, reset gene data
                gene_data = [line]
                current_gene = gene

    if current_gene in genes:
        yield gene_data



def create_annotation_file(annot,out_annot):
    """
    Creates regenie-ready annotation file
    """
    with open(annot,'rt') as i, open(out_annot,'wt') as o :
        line = next(i)
        for line in i:
            chrom,pos,ref,alt,gene_most_severe,most_severe = line.strip().split()
            o.write('\t'.join([f"chr{chrom}_{pos}_{ref}_{alt}",gene_most_severe,most_severe]) + '\n')


def create_sets(annot,set_file):
    """
    Builds sets for regenie

    """
    gene_dict = defaultdict(list)

    with open(annot) as i:
        for line in i:
            snp,gene,most_severe = line.strip().split()
            gene_dict[gene].append(snp)
    print('gene dict built')
    with open(set_file,'wt') as o:
        for gene in gene_dict:
            variants = gene_dict[gene]
            chrom_list,pos_list = zip(*[elem.split("_")[:2] for elem in variants])
            chrom_set = set(chrom_list)
            if len(chrom_set) ==1:
                chrom = chrom_list[0].split("chr")[1]
                o.write('\t'.join([gene,chrom,pos_list[0],','.join(variants)]) +'\n')
            else:
                print(gene,chrom_set)


def check_error_regenie(stdout):
    return "ERROR" in stdout or "End time" not in stdout


def read_phenos(pred):
    with open(pred) as i: phenos = [elem.strip().split()[0] for elem in i]
    return phenos

def get_param_list(out,name,pred,chrom_list,phenos,force):
    """
    Get list of missing regenie logs
    """

    param_list = list(product(phenos,chrom_list))

    if force:
        return param_list

    else:
        params = []
        for pheno,chrom in param_list:
            log_file = os.path.join(out,f"{name}_{chrom}_{pheno}.log")
            regenie_file = os.path.join(out,f"{name}_{chrom}_{pheno}.regenie")
            check_run = True

            # if force it needs to be run
            if force:
                check_run = False
            # if the regenie output file doesn't exist or is empty
            if os.path.isfile(regenie_file):
                if not mapcount(regenie_file):
                    check_run = False
            else:
                check_run = False
            # if the log exists and it doesn't show signs of failure
            if os.path.isfile(log_file):
                with open(log_file) as i: stdout = i.read()
                if check_error_regenie(stdout):
                    check_run = False
            else:
                check_run = False

            if not check_run:
                params.append((chrom,pheno))
        print(f"Run params: {params}")
        return params
