import os
import json
import numpy as np
from pyjaspar import jaspardb
from esdeg.functions import get_threshold, calculate_gc, make_k_mers, jaspar_to_pwm
from esdeg.parsers import promoters_parser
import esdeg.speedup as sup


def prepare_motif_db(output_dir, path_to_promoters, taxon):

    #create outdir
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    print('-'*30)
    print('Read promoters')
    order = 0
    promoters, ids = promoters_parser(path_to_promoters)
    gc_content = calculate_gc(promoters)
    kmer_converter = make_k_mers(order)
    promoters = sup.seq_to_int(promoters, order, kmer_converter)
    print('-'*30)

    #metadata
    metadata = dict()
    metadata['ids'] = list(ids)
    metadata['gc'] = list(gc_content)
    metadata['motif_id'] = []
    metadata['tf_name'] = []
    metadata['tf_class'] = []

    print('Read matrices')
    #plants
    #vertebrates
    #insects
    #urochordates
    #nematodes
    #fungi
    jdb_obj = jaspardb(release='JASPAR2022')
    motifs = jdb_obj.fetch_motifs(collection = 'CORE',tax_group = [taxon], min_length=8)
    number_of_motifs = len(motifs)
    print(f'Number of matrices = {number_of_motifs}')
    print('-'*30)

    for index, motif_data in enumerate(motifs, start=1):
        motif_id = motif_data.matrix_id
        tf_name = motif_data.name
        tf_class = motif_data.tf_class
        if len(tf_class) == 0:
            tf_class = 'NA'
        elif len(tf_class) == 1:
            tf_class = tf_class[0]
        else:
            tf_class = '::'.join(tf_class)

        pwm = jaspar_to_pwm(motif_data)

        metadata['motif_id'].append(motif_id)
        metadata['tf_name'].append(tf_name)
        metadata['tf_class'].append(tf_class)

        print(f'{index}. {motif_id} - {tf_name}')
        scores = sup.scaner_pwm(promoters, pwm)
        flatten_scores = scores.ravel()
        flatten_scores = sup.sort(flatten_scores)
        threshold_table = get_threshold(flatten_scores)
        threshold_table = np.array(threshold_table)
        fprs_table = threshold_table[:,1]
        fprs_choosen = np.power(10, np.linspace(np.log10(5e-4), np.log10(1e-4), num=30))
        indexes = np.searchsorted(fprs_table, fprs_choosen)
        threshold_table = threshold_table[indexes]
        counts = sup.count_sites(scores, threshold_table)
        np.save(f'{output_dir}/{motif_id}.npy', counts)
    with open(f'{output_dir}/metadata.json', 'w') as file:
        json.dump(metadata, file)
    pass
