import os
import numpy as np
import gzip
from collections import defaultdict as dd
import pickle
import scipy
import shutil
import time
import sys
import mmap
import multiprocessing
import subprocess
from subprocess import call
import shlex
import csv
cpus = multiprocessing.cpu_count()

rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
# REQUIRED FILES

def get_progress_test(args):
    

    import glob
    matrix_chunk_path = os.path.join(args.out_path,'matrix_chunks/')
    chunk_path = os.path.join(matrix_chunk_path, '*gene_matrix*')
    files = glob.iglob(chunk_path)
    lines = 0
    for f in files:
        if os.path.getsize(f) > 0:
            lines += mapcount(f)
    
    progressBar(lines,len(args.genes))

def make_sure_path_exists(path):
    import errno
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise                


tmp_path = rootPath + 'tmp/'
make_sure_path_exists(tmp_path)
plink = shutil.which('plink')
plink2 = shutil.which('plink2')



mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
mem_mib = mem_bytes/(1024.**2) 
proc_mem = mem_mib / (cpus +1)




def split_array_chunk(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def get_progress(args):
    
    if args.test:
        genes = args.test * args.cpus
    else:
        with open(args.g2v,'rb') as i: g2v = pickle.load(i)
        gene_list = np.array(list(g2v.keys()))
        genes = len(gene_list)

    import glob    
    chunk_path = os.path.join(args.tmp_path, args.lof +'_' + 'matrix_chunk*')
    files = glob.iglob(chunk_path)
    lines = 0
    for f in files:
        if os.path.getsize(f) > 0:
            lines += mapcount(f)
    
    progressBar(lines,genes)

def get_path_info(path):
    file_path = os.path.dirname(path)
    basename = os.path.basename(path)
    file_root, file_extension = os.path.splitext(basename)
    return file_path,file_root,file_extension


def identify_separator(line):
    import csv
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(line)
    return dialect.delimiter

def return_header_variants(matrixPath):
    '''
    The plink command adds the alt to the name of the variant. Here i loop through the variants and just return the original variant name
    '''
    with open(matrixPath,'rt') as i:header = i.readline().strip()
    sep = identify_separator(header)
    header = header.split(sep)
    header_variants = ['_'.join(elem.split('_')[:-1]) for elem in header]
        
    return header_variants

def pad(s,padding = ' '):
    return padding +s + padding


def tmp_bash(cmd):
    from tempfile import NamedTemporaryFile

    scriptFile = NamedTemporaryFile(delete=True)
    with open(scriptFile.name, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(cmd + "\n")

    os.chmod(scriptFile.name,0o777)
    scriptFile.file.close()
    subprocess.check_call(scriptFile.name)

def run_bash_script(cmd,file_path):
    with open(file_path,'wt') as o:
        o.write('#!/bin/bash\n')
        o.write(cmd)
    call(shlex.split('chmod +x ' + file_path))
    call(file_path,shell = True)
    
def return_open_func(f):
    '''
    Detects file extension and return proper open_func
    '''
    import gzip
    from functools import partial

    file_path,file_root,file_extension = get_path_info(f)
    print(file_extension)

    if 'bgz' in file_extension:
        print('gzip.open with rb mode')
        open_func = partial(gzip.open, mode = 'rb')
    
    elif 'gz' in file_extension:
        print('gzip.open with rt mode')
        open_func = partial(gzip.open, mode = 'rt')

    else:
        print('regular open')
        open_func = open      
    return open_func

def progressBar(value, endvalue, bar_length=20):
    '''
    Writes progress bar, ven value (eg.current row) and endvalue(eg. total number of rows)
    '''

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()
    
def timing_function(some_function):

    """
    Outputs the time a function takes  to execute.
    """

    def wrapper(*args,**kwargs):
        t1 = time.time()
        some_function(*args)
        t2 = time.time()
        print("Time it took to run the function: " + str((t2 - t1)))

    return wrapper



def file_exists(fname):
    '''
    Function to pass to type in argparse
    '''
    if os.path.isfile(fname):
        return str(fname)
    else:
        print('File does not exist')
        sys.exit(1)



    

def pretty_print(string,l = 30):
    l = l-int(len(string)/2)
    print('-'*l + '> ' + string + ' <' + '-'*l)
    


def mapcount(filename):
    '''
    Counts line in file
    '''
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines

class SliceMaker(object):
    '''
    allows to pass around slices
    '''
    def __getitem__(self, item):
        return item

def line_iterator(f,separator ='\t',count = False,columns = SliceMaker()[:],dtype = str,skiprows = 0):
    '''
    Function that iterates through a file and returns each line as a list with separator being used to split.
    N.B. it requires that all elements are the same type
    '''
    if f.split('.')[-1] != 'gz':
        i = open(f,'rt')

    elif f.split('.')[-1] == 'gz':
        i = gzip.open(f,'rt')

    for x in range(skiprows):next(i)

    if count is False:
        for line in i:
            yield np.array(line.strip().split(separator),dtype = str)[columns].astype(dtype)
    else:
        row = 0
        for line in i:
            row +=1 
            yield row,np.array(line.strip().split(separator),dtype = str)[columns].astype(dtype)



def basic_iterator(f,separator ='\t',skiprows = 0,count = False,columns = 'all'):
    '''
    Function that iterates through a file and returns each line as a list with separator being used to split.
    '''
    if f.split('.')[-1] != 'gz':
        i = open(f,'rt')

    elif f.split('.')[-1] == 'gz':
        i = gzip.open(f,'rt')

    for x in range(skiprows):next(i)

    if count is False:
        for line in i:
            line =line.strip().split(separator)
            line = return_columns(line,columns)
            yield line
    else:
        row = 0
        for line in i:
            line =line.strip().split(separator)
            line = return_columns(line,columns)
            row += 1   
            yield row,line
def return_columns(l,columns):
    '''
    Returns all columns, or rather the elements, provided the columns
    '''
    if columns == 'all':
        return l
    elif type(columns) == int:
        return l[columns]
    elif type(columns) == list:
        return list(map(l.__getitem__,columns))
