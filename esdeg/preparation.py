import os
import json
import numpy as np
from esdeg.functions import get_threshold, calculate_gc, make_k_mers
from esdeg.parsers import matrix_parser, bamm_parser, promoters_parser
import esdeg.speedup as sup


def prepare_motif_db_pwm(path_to_db, output_dir, path_to_promoters, organism, file_format='meme'):
    
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
    metadata['organism'] = organism
    metadata['ids'] = list(ids)
    metadata['gc'] = list(gc_content)
    metadata['matrices'] = []
    
    print('Read matrices')
    motifs = matrix_parser(path_to_db, f=file_format)
    number_of_motifs = len(motifs)
    print(f'Number of matrices = {number_of_motifs}')
    print('-'*30)
    
    for index, motif_data in enumerate(motifs, start=1):
        name, pwm, pfm, motif_length = motif_data
        metadata['matrices'].append(name)
        print(f'{index}. {name}')
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
        np.save(f'{output_dir}/{name}.npy', counts)
    with open(f'{output_dir}/metadata.json', 'w') as file:
        json.dump(metadata, file)
    pass


def prepare_motif_db_bamm(path_to_db, order, output_dir, path_to_promoters, organism):
    
    #create outdir
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    print('-'*30)
    print('Read promoters')
    promoters, ids = promoters_parser(path_to_promoters)
    gc_content = calculate_gc(promoters)
    kmer_converter = make_k_mers(order)
    promoters = sup.seq_to_int(promoters, order, kmer_converter)
    print('-'*30)

    
    #metadata
    metadata = dict()
    metadata['organism'] = organism
    metadata['ids'] = list(ids)
    metadata['gc'] = list(gc_content)
    metadata['matrices'] = []
    
    print('Read motifs (BaMM)')
    motifs = bamm_parser(path_to_db)
    number_of_motifs = len(motifs)
    print(f'Number of matrices = {number_of_motifs}')
    print('-'*30)
    
    for index, motif_data in enumerate(motifs, start=1):
        name, bamm, m_order, motif_length = motif_data
        metadata['matrices'].append(name)
        print(f'{index}. {name}')
        scores = sup.scaner_bamm(promoters, bamm, order)
        flatten_scores = scores.ravel()
        flatten_scores = sup.sort(flatten_scores)
        threshold_table = get_threshold(flatten_scores)
        threshold_table = np.array(threshold_table)
        fprs_table = threshold_table[:,1]
        fprs_choosen = np.power(10, np.linspace(np.log10(5e-4), np.log10(1e-4), num=30))
        indexes = np.searchsorted(fprs_table, fprs_choosen)
        threshold_table = threshold_table[indexes]
        counts = sup.count_sites(scores, threshold_table)
        np.save(f'{output_dir}/{name}.npy', counts)
    with open(f'{output_dir}/metadata.json', 'w') as file:
        json.dump(metadata, file)
    pass