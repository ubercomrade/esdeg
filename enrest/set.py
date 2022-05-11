import pandas as pd
import numpy as np
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from functools import partial
from enrest.functions import *
from enrest.parsers import matrices_parser, promoters_parser, read_set_of_genes
from enrest.scaner import scaner
import enrest.speedup as sup


def work_with_matrix(args, set_ids=None, all_ids=None, promoters=None, parameter=None):
    name, pwm, pfm, matrix_length, middle_score = args
    line = {('', 'ID'): name}
    print(f'{name}')
    all_scores = scaner(promoters, pwm)
    best_scores = np.max(all_scores)
    flatten_scores = all_scores.ravel()
    flatten_scores = sup.sort(flatten_scores)
    threshold_table = get_threshold(flatten_scores)
    threshold_table = np.array(threshold_table)
    fprs_table = threshold_table[:,1]
    fprs_choosen = np.array([0.0005, 0.00015, 0.00005]) # LOW, MIDDLE, HIGH
    indexes = np.searchsorted(fprs_table, fprs_choosen)
    threshold_table = threshold_table[indexes]
    other_ids = get_other_gene_ids_for_set_case(set_ids, all_ids)
    if parameter == "enrichment":
        set_scores, other_scores, genes = split_scores_by_gene_ids(all_scores, all_ids, set_ids, other_ids)
    elif parameter == "fraction":
        set_scores, other_scores, genes = split_scores_by_gene_ids(best_scores, all_ids, set_ids, other_ids)
    results = run_test(genes, set_scores, other_scores, threshold_table, parameter)
    line.update(results)
    return line


def set_case(path_to_set, path_to_db, output_dir, path_to_promoters, 
             file_format="meme", parameter="enrichment", number_of_cores=2):    
    print('-'*30)
    print('Read SET of genes')
    set_ids = read_set_of_genes(path_to_set)
    print('-'*30)
    print('Read promoters')
    promoters, all_ids = promoters_parser(path_to_promoters)
    print('-'*30)
    print('Read matrices')
    matrices = matrices_parser(path_to_db, f=file_format)
    number_of_matrices = len(matrices)
    print(f'Number of matrices = {number_of_matrices}')
    print('-'*30)
    with Pool(number_of_cores) as pool:
        results = pool.map(partial(work_with_matrix, set_ids=set_ids, all_ids=all_ids, promoters=promoters, parameter=parameter), matrices)
        results = list(results)
    df = pd.DataFrame(results, columns=results[0].keys())
    output_path = f"{output_dir}/all.tsv"
    df.to_csv(output_path, sep='\t', index=False)
    print('-'*30)
    print('All done. Exit')
    return None
