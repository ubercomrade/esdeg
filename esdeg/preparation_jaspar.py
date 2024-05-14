import os
import json
import numpy as np
from itertools import repeat
from collections import Counter
from pyjaspar import jaspardb
from esdeg.functions import get_threshold, calculate_gc, \
jaspar_to_pwm, seq_to_int, all_scores_pwm, count_sites
from esdeg.parsers import promoters_parser
from concurrent.futures import ProcessPoolExecutor


def check_motif_annotaion(ann):
    if len(ann) == 0:
        ann = 'NA'
    elif len(ann) == 1:
        ann = ann[0]
    else:
        ann = '::'.join(ann)
    return ann

def worker(index, motif_data, promoters, output_dir):

    motif_id = motif_data.matrix_id
    tf_name = motif_data.name
    tf_class = check_motif_annotaion(motif_data.tf_class)
    tf_family = check_motif_annotaion(motif_data.tf_family)

    pwm = jaspar_to_pwm(motif_data)
    pwm = np.concatenate((pwm, np.array([pwm.shape[1] * [-100]])), axis=0)

    print(f'{index}. {motif_id} - {tf_name}')
    scores = all_scores_pwm(promoters, pwm)
    flatten_scores = scores.ravel()
    flatten_scores = np.sort(flatten_scores)[::-1]
    threshold_table = get_threshold(flatten_scores)
    threshold_table = np.array(threshold_table)
    fprs_table = threshold_table[:,1]
    fprs_choosen = np.power(10, np.linspace(np.log10(5e-5), np.log10(1e-3), num=20))
    indexes = np.searchsorted(fprs_table, fprs_choosen)
    threshold_table = threshold_table[indexes]
    actual_fprs = threshold_table[:,1]
    number_of_uniq_fprs = len(np.unique(actual_fprs))
    counts = count_sites(scores, threshold_table)
    np.save(f'{output_dir}/{motif_id}.npy', counts)

    motif_info = {
        'motif_id': motif_id,
        'tf_name': tf_name,
        'tf_class': tf_class,
        'tf_family': tf_family,
        'number_of_uniq_fprs': number_of_uniq_fprs
    }
    return motif_info


def prepare_motif_db(output_dir, path_to_promoters, taxon, nproc):

    #create outdir
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    print('-'*30)
    print('Read promoters')
    promoters, ids = promoters_parser(path_to_promoters)

    # Filter promoters if they have different length
    length_of_promoters = [len(i) for i in promoters]
    length_dict = Counter(length_of_promoters)
    if len(length_dict) > 1:
        print('#'*10 + '  WARNING!  ' + '#'*10)
        print(f'Fasta file contain sequnces with different legth')
        for index, (length, counts) in enumerate(length_dict.items(), 1):
            print(f'{index}. {counts} number of sequnces with length {length}')
        length, counts = length_dict.most_common(1)[0]
        print(f'Only sequnces with most common length (={length}) will be used')
        filter = [True if i == length else False for i in length_of_promoters]
        promoters = [i for i, logical in zip(promoters, filter) if logical]
        ids = np.array([i for i, logical in zip(ids, filter) if logical])
        print('#'*10 + '  DATA FILTERED!  ' + '#'*10)
    # End filter

    gc_content = calculate_gc(promoters)
    promoters = seq_to_int(promoters)
    print('-'*30)

    #metadata
    metadata = dict()
    metadata['db'] = 'jaspar'
    metadata['taxon'] = taxon
    metadata['ids'] = list(ids)
    metadata['gc'] = list(gc_content)
    metadata['motif_id'] = []
    metadata['tf_name'] = []
    metadata['tf_class'] = []
    metadata['tf_family'] = []

    print('Read matrices')
    #plants
    #vertebrates
    #insects
    #urochordates
    #nematodes
    #fungi
    jdb_obj = jaspardb(release='JASPAR2024')
    motifs = jdb_obj.fetch_motifs(collection = 'CORE',tax_group = [taxon], min_length=6)
    number_of_motifs = len(motifs)
    print(f'Number of matrices = {number_of_motifs}')
    print('-'*30)

    with ProcessPoolExecutor(max_workers=nproc) as executor:
        information_of_motifs = executor.map(worker, range(1, number_of_motifs + 1), motifs, repeat(promoters), repeat(output_dir))

    for motif_info in information_of_motifs:
        for i in ['motif_id', 'tf_name', 'tf_class', 'tf_family']:
            metadata[i].append(motif_info[i])
            if i == 'motif_id':
                metadata[motif_info[i]] = {'number_of_uniq_fprs': motif_info['number_of_uniq_fprs']}

    with open(f'{output_dir}/metadata.json', 'w') as file:
        json.dump(metadata, file)

    pass
