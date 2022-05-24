import pandas as pd
import numpy as np
from enrest.functions import run_test, get_threshold
from enrest.parsers import matrices_parser, promoters_parser, read_set_of_genes
import enrest.speedup as sup


def work_with_matrix(name, pwm, pfm, matrix_length, all_ids, deg_table, promoters, parameter, 
    padj_thr, log2fc_thr_deg, log2fc_thr_background):
    print(f'{name}')
    container = {'ALL': [], 'UP': [], 'DOWN': []}
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
    for index, condition in enumerate(['ALL', 'UP', 'DOWN'], 1):
        line = {('', 'ID'): name}
        deg_ids = get_deg_gene_ids(deg_table, condition, padj_thr=padj_thr, log2fc_thr=log2fc_thr_deg)
        other_ids = get_other_gene_ids_for_deg_case(deg_table, padj_thr=padj_thr, log2fc_thr=log2fc_thr_background)
        if parameter == "enrichment":
            deg_scores, other_scores, genes = split_scores_by_gene_ids(all_scores, all_ids, deg_ids, other_ids)
        elif parameter == "fraction":
            deg_scores, other_scores, genes = split_scores_by_gene_ids(best_scores, all_ids, deg_ids, other_ids)
        results = run_test(genes, deg_scores, other_scores, threshold_table, parameter)
        line.update(results)
        container[condition] = line
    return container


def deg_case(path_to_deg, path_to_db, output_dir, path_to_promoters, 
             file_format='meme', parameter='enrichment', padj_thr=0.05,
             log2fc_thr_deg=1, log2fc_thr_background=np.log2(5/4), number_of_cores=2):    
    print('-'*30)
    print('Read DEG table')
    deg_table = pd.read_csv(path_to_deg, sep=',', comment='#')
    deg_table = deg_table[deg_table['padj'] <= 1]
    print('-'*30)
    print('Read promoters')
    promoters, all_ids = promoters_parser(path_to_promoters)
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
                         all_ids, deg_table, 
                         promoters, parameter, padj_thr, 
                         log2fc_thr_deg, log2fc_thr_background)
        results.append(line)
    for index, condition in enumerate(['ALL', 'UP', 'DOWN'], 1):
        container = [i[condition] for i in results]
        df = pd.DataFrame(container, columns=container[0].keys())
        condition = condition.lower()
        output_path = f"{output_dir}/{condition}.tsv"
        df.to_csv(output_path, sep='\t', index=False)
    print('-'*30)
    print('All done. Exit')
    return None
