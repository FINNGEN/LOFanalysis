import os.path
from Tools.utils import mapcount,progressBar

def get_progress_test(args):
    

    import glob
    
    matrix_chunk_path = os.path.join(args.out_path,'matrix_chunks/',args.lof + '_CHUNK_gene_matrix.tsv')

    files = [matrix_chunk_path.replace('CHUNK',str(i)) for i in range(args.cpus)]
    lines = 0
    for f in files:
        if os.path.getsize(f) > 0:
            lines += mapcount(f)
    
    progressBar(lines,len(args.genes))



def split_array_chunk(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out
