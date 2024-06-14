import sys
import json
import pandas as pd
import numpy as np
from importlib import resources
from esdeg.functions import run_test, get_deg_gene_ids, get_other_gene_ids_for_deg_case, split_by_gene_ids, get_motif_to_cluster, fdrcorrection_log


def deg_case(path_to_deg, path_to_db,
             parameter='enrichment',
             padj_thr=0.05,
             log2fc_thr_deg=1,
             log2fc_thr_background=np.log2(5/4),
             gc_threshold=0.25,
             condition='down'):

    print('Read metadata')
    file = open(f'{path_to_db}/metadata.json')
    metadata = json.load(file)
    file.close()
    gc_content = np.array(metadata['gc'])
    ids = np.array(metadata['ids'])
    taxon = metadata['taxon']
    print('-'*30)

    print('Read DEG table')
    deg_table = pd.read_csv(path_to_deg, sep=',', comment='#')
    deg_table = deg_table[deg_table['padj'] <= 1]
    foreground_ids = get_deg_gene_ids(deg_table, condition, padj_thr=padj_thr, log2fc_thr=log2fc_thr_deg)
    if len(foreground_ids) <= 5:
        ll = len(foreground_ids)
        print(f'WARNING: Namber of DEGs ({ll}) is dramaticaly low!')
    check_intersection = np.sum(np.in1d(foreground_ids, ids, assume_unique=True))
    if check_intersection == 0:
        print(f'There are no common IDs in the database with the DEG IDs. Maybe different identifiers are used. Exit.')
        sys.exit()
    other_ids = get_other_gene_ids_for_deg_case(deg_table, padj_thr=padj_thr, log2fc_thr=log2fc_thr_background)
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
