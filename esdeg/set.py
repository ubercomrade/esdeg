import sys
import json
import pandas as pd
import numpy as np
from importlib import resources
from esdeg.functions import run_test, get_other_gene_ids_for_set_case, split_by_gene_ids, get_motif_to_cluster, fdrcorrection_log
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
    taxon = metadata['taxon']
    print('-'*30)

    print('Read SET of genes')
    foreground_ids = read_set_of_genes(path_to_set)
    check_intersection = np.sum(np.in1d(foreground_ids, ids, assume_unique=True))
    if check_intersection == 0:
        print(f'There are no common IDs in the database with the DEG IDs. Maybe different identifiers are used. Exit.')
        sys.exit()
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
        tf_family = metadata['tf_family'][index - 1]
        number_of_uniq_fprs = metadata[motif_id]['number_of_uniq_fprs']
        line = {'motif_id': motif_id,
                'tf_name': tf_name,
                'tf_class': tf_class,
                'tf_family': tf_family}
        print(f'{index} {motif_id} - {tf_name}')
        with np.load(f'{path_to_db}/{motif_id}.npz') as file:
            counts = file['counts']
        foreground, foreground_gc, other, other_gc, genes = split_by_gene_ids(counts,
                                                                              gc_content,
                                                                              ids,
                                                                              foreground_ids,
                                                                              other_ids)
        out, log_pvalues, odds_ratios, index_of_best = run_test(genes, foreground,
                                                          foreground_gc, other,
                                                          other_gc, gc_threshold,
                                                          parameter, number_of_uniq_fprs)
        line.update(out)
        results.append(line)
    df = pd.DataFrame(results)
    lg_adj_pval = fdrcorrection_log(df['ln(pval)']) # natrual log
    df['adj.pval'] = np.exp(lg_adj_pval)
    df['log10(pval)'] = df['ln(pval)'] / np.log(10) # from loge to log10
    df['log10(adj.pval)'] = lg_adj_pval / np.log(10) # from loge to log10

    if metadata['db'] == 'jaspar':
        cluster_path = resources.files('esdeg').joinpath(f'clusters/{taxon}.tsv')
        motif_to_cluster = get_motif_to_cluster(cluster_path)
        df['jaspar_cluster'] = [motif_to_cluster[i] for i in df['motif_id']]

    #split dimers (JASPAR)
    if metadata['db'] == 'jaspar':
        df['tf_name'] = df['tf_name'].str.split('::')
        df['tf_class'] = df['tf_class'].str.split('::')
        df['tf_family'] = df['tf_family'].str.split('::')


        # fix number of Class and Family (reason is MA1628.2, JASPAR bug)
        container = []
        for index, i in df.iterrows():
            if len(i['tf_name']) != len(i['tf_class']) and len(i['tf_class']) != len(i['tf_family']):
                length = np.max([len(i['tf_name']), len(i['tf_class']), len(i['tf_family'])])
                if len(i['tf_name']) != length:
                    i['tf_name'] = length * i['tf_name']
                if len(i['tf_class']) != length:
                    i['tf_class'] = length * i['tf_class']
                if len(i['tf_family']) != length:
                    i['tf_family'] = length * i['tf_family']
            container.append(i)
        df = pd.DataFrame(container)
        # end fix

        df = df.explode(column=['tf_name', 'tf_class', 'tf_family'], ignore_index=True)
        df = df.sort_values(by='adj.pval')
    # end split

    print('-'*30)
    return df, taxon
