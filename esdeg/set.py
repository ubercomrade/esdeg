import json
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from esdeg.functions import run_test, get_other_gene_ids_for_set_case, split_by_gene_ids
from esdeg.parsers import read_set_of_genes


def set_case(path_to_set, path_to_db,
             parameter="enrichment",
             gc_threshold=0.25):

    print('Read metadata')
    file = open(f'{path_to_db}/metadata.json')
    metadata = json.load(file)
    file.close()
    gc_content = np.array(metadata['gc'])
    ids = np.array(metadata['ids'])

    print('-'*30)
    print('Read SET of genes')
    foreground_ids = read_set_of_genes(path_to_set)
    check_intersection = np.sum(np.in1d(foreground_ids, ids, assume_unique=True))
    if check_intersection == 0:
        print(f'There are no common IDs in the database with the DEG IDs. Maybe different identifiers are used. Exit.')
        os.exit()
    other_ids = get_other_gene_ids_for_set_case(foreground_ids, ids)
    print('-'*30)


    print('Work with matrices')
    number_of_matrices = len(metadata['motif_id'])
    print(f'Number of matrices = {number_of_matrices}')
    print('-'*30)

    results = []
    for index, motif_id in enumerate(metadata['motif_id'], 1):
        tf_name = metadata['tf_name'][index - 1]
        tf_class = metadata['tf_class'][index - 1]
        line = {'motif_id': motif_id,
                'tf_name': tf_name,
                'tf_class': tf_class}
        print(f'{index} {motif_id} - {tf_name}')
        counts = np.load(f'{path_to_db}/{motif_id}.npy')
        foreground, foreground_gc, other, other_gc, genes = split_by_gene_ids(counts,
                                                                              gc_content,
                                                                              ids,
                                                                              foreground_ids,
                                                                              other_ids)
        out = run_test(genes, foreground, foreground_gc, other, other_gc, gc_threshold, parameter)
        line.update(out)
        results.append(line)
    df = pd.DataFrame(results)
    _, adj_pval, _, _ = multipletests(df['pval'], method='fdr_bh')
    df['adj.pval'] = adj_pval
    print('-'*30)
    return df
