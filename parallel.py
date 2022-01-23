from joblib import Parallel, delayed
import multiprocessing
from itertools import islice
import os
import math


def chunks(data, num_seqs):
    # Create batches
    it = iter(data)
    for i in range(0, len(data), num_seqs):
        yield {k: data[k] for k in islice(it, num_seqs)}


def rnafold(file_name):
    cmd = f'RNAfold -p -d2 --noLP --noPS --noDP < {file_name} > {file_name[:-3]}_fold.out'
    os.system(cmd)


def rnafold_in_parallel(iso2, output_file_stem, ac):
    num_cores = multiprocessing.cpu_count() - 1
    chunk_size = math.ceil(len(iso2.trnas)/num_cores)
    trna_batches = chunks(iso2.trnas, chunk_size)

    file_names = []
    for i, batch in enumerate(trna_batches, 1):
        with open(f'{output_file_stem}_{ac}_{i}.fa', 'w+') as f:
            for name, trna in batch.items():
                f.write(f'>{name}\n{trna.seq[ac]}\n')
        file_names.append(f'{output_file_stem}_{ac}_{i}.fa')

    Parallel(n_jobs=num_cores)(delayed(rnafold)(filename) for filename in file_names)

    outfile_names = [f'{file_name[:-3]}_fold.out' for file_name in file_names]
    with open(f'{output_file_stem}_{ac}_complete_fold.out', 'w') as outfile:
        for f_name in outfile_names:
            with open(f_name) as infile:
                for line in infile:
                    outfile.write(line)
