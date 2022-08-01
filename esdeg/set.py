import json
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from esdeg.functions import run_test, get_other_gene_ids_for_set_case, split_by_gene_ids
from esdeg.parsers import read_set_of_genes


def set_case(path_to_set, path_to_db, path_to_output, organism, 
             parameter="enrichment", gc_threshold=0.25): 
    
    print('Read metadata')
    file = open(f'{path_to_db}/metadata.json')
    metadata = json.load(file)
    file.close()
    if metadata['organism'] != organism:
        md_org = metadata['organism']
        print(f'DB was prepared for {md_org}, but you try to use it for {organism}. Exit.')
        os.exit()
    gc_content = np.array(metadata['gc'])
    ids = np.array(metadata['ids'])

    print('-'*30)
    print('Read SET of genes')
    foreground_ids = read_set_of_genes(path_to_set)
    other_ids = get_other_gene_ids_for_set_case(foreground_ids, ids)
    print('-'*30)
    
    
    print('Work with matrices')
    number_of_matrices = len(metadata['matrices'])
    print(f'Number of matrices = {number_of_matrices}')
    print('-'*30)

    results = []
    for index, matrix_name in enumerate(metadata['matrices'], 1):
        line = {'motif': matrix_name}
        print(f'{index} {matrix_name}')
        counts = np.load(f'{path_to_db}/{matrix_name}.npy')
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
    df = df[['motif', 'log(or)', 'distance', 'pval', 'adj.pval', 'genes']]
    df.to_csv(path_to_output, sep='\t', index=False)
    print('-'*30)
    print('All done. Exit')
    return None
