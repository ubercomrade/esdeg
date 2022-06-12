import pandas as pd
import numpy as np
from enrest.functions import run_test_fasta, get_threshold, calculate_gc
from enrest.parsers import matrices_parser, fasta_parser, promoters_parser, read_set_of_genes
import enrest.speedup as sup


def work_with_matrix(name, pwm, pfm, matrix_length, 
                     foreground, foreground_gc,
                     background, background_gc, 
                     promoters, parameter, gc_threshold):
    line = {('', 'ID'): name}
    print(f'{name}')
    all_scores = sup.scaner(promoters, pwm)
    best_scores = np.max(all_scores, axis=1)
    flatten_scores = all_scores.ravel()
    flatten_scores = sup.sort(flatten_scores)
    threshold_table = get_threshold(flatten_scores)
    threshold_table = np.array(threshold_table)
    fprs_table = threshold_table[:,1]
    fprs_choosen = np.array([0.0005, 0.00015, 0.00005]) # LOW, MIDDLE, HIGH
    indexes = np.searchsorted(fprs_table, fprs_choosen)
    threshold_table = threshold_table[indexes]
    if parameter == "enrichment":
        foreground_scores = sup.scaner(foreground, pwm)
        background_scores = sup.scaner(background, pwm)
    elif parameter == "fraction":
        foreground_scores = np.max(sup.scaner(foreground, pwm), axis=1)
        background_scores = np.max(sup.scaner(background, pwm), axis=1)
    results = run_test_fasta(foreground_scores, foreground_gc,
                             background_scores, background_gc,
                             threshold_table, gc_threshold, parameter)
    line.update(results)
    return line


def fasta_case(path_to_foreground, path_to_background, path_to_db, output_dir, path_to_promoters, 
             file_format="meme", parameter="enrichment", gc_threshold=0.1):    
    print('-'*30)
    print('Read fasta foreground and background')
    foreground = fasta_parser(path_to_foreground)
    foreground = sup.seq_to_int(foreground)
    foreground_gc = calculate_gc(foreground)
    
    background = fasta_parser(path_to_background)
    background = sup.seq_to_int(background)
    background_gc = calculate_gc(background)
    print('-'*30)
    print('Read promoters')
    promoters, all_ids = promoters_parser(path_to_promoters)
    promoters = sup.seq_to_int(promoters)
    print('-'*30)
    print('Read matrices')
    matrices = matrices_parser(path_to_db, f=file_format)
    number_of_matrices = len(matrices)
    print(f'Number of matrices = {number_of_matrices}')
    print('-'*30)
    results = []
    for matrix_data in matrices:
        name, pwm, pfm, matrix_length = matrix_data
        line = work_with_matrix(name, pwm, pfm, matrix_length, 
                     foreground, foreground_gc,
                     background, background_gc, 
                     promoters, parameter, gc_threshold)
        results.append(line)
    df = pd.DataFrame(results, columns=results[0].keys())
    output_path = f"{output_dir}/all.tsv"
    df.to_csv(output_path, sep='\t', index=False)
    print('-'*30)
    print('All done. Exit')
    return None
