import pandas as pd
import numpy as np


def create_gname_converter(gtf_path):
    gtf = pd.read_csv(gtf_path, sep='\t', comment='#', dtype = {0: str}, header=None)
    gtf = gtf[gtf[2] == 'gene']
    gtf = gtf[[True if "gene_name" in i else False for i in gtf[8]]]
    converter = []
    for line in gtf.itertuples():
        gene_id = line[9].split(';')[0].split()[1].strip('\"')
        gene_name = line[9].split(';')[2].split()[1].strip('\"')
        converter.append((gene_id, gene_name.upper()))
    #name_to_id = {i[1]:i[0] for i in converter}
    id_to_name = {i[0]:i[1] for i in converter}
    return id_to_name


def read_esdeg_table(esdeg_path):
    esdeg = pd.read_csv(esdeg_path, sep='\t')
    esdeg = esdeg[esdeg['tf_family'] != 'NA']
    esdeg['tf_name'] = esdeg['tf_name'].str.upper()
    esdeg = esdeg.reset_index(drop=True)

    # #split dimers (JASPAR)
    # esdeg['tf_name'] = esdeg['tf_name'].str.split('::')
    # esdeg['tf_class'] = esdeg['tf_class'].str.split('::')
    # esdeg['tf_family'] = esdeg['tf_family'].str.split('::')
    #
    # # fix number of Class and Family (reason is MA1628.2, JASPAR bug)
    # container = []
    # for index, i in esdeg.iterrows():
    #     if len(i['tf_name']) != len(i['tf_class']) and len(i['tf_class']) != len(i['tf_family']):
    #         length = np.max([len(i['tf_name']), len(i['tf_class']), len(i['tf_family'])])
    #         if len(i['tf_name']) != length:
    #             i['tf_name'] = length * i['tf_name']
    #         if len(i['tf_class']) != length:
    #             i['tf_class'] = length * i['tf_class']
    #         if len(i['tf_family']) != length:
    #             i['tf_family'] = length * i['tf_family']
    #     container.append(i)
    # esdeg = pd.DataFrame(container)
    # # end fix
    #
    # esdeg = esdeg.explode(column=['tf_name', 'tf_class', 'tf_family'], ignore_index=True)
    # esdeg = esdeg.sort_values(by='adj.pval')
    # # end split

    #esdeg = esdeg.drop_duplicates(subset=['tf_name'], keep='first')
    #esdeg = esdeg.set_index(['tf_name'])
    esdeg = esdeg.rename(columns={"adj.pval": "me_padj"})
    esdeg = esdeg.drop(columns=['log10(pval)', 'log10(adj.pval)'])
    return esdeg


def annotation(deg_path, counts_path, esdeg_path, gtf_path, filter_flag=False, best_flag=False,
               me_padj_thr=0.05, de_padj_thr=0.05, lor_thr=1., lfc_thr=1., counts_filter=5):


    id_to_name = create_gname_converter(gtf_path)
    esdeg = read_esdeg_table(esdeg_path)

    # DEGs
    deg = pd.read_csv(deg_path, sep=',')
    deg = deg[[True if i in id_to_name else False for i in deg['id']]]
    deg = deg.assign(tf_name = [id_to_name[i] for i in deg['id']])
    deg = deg.set_index(['tf_name'])
    # DEGs


    # counts
    counts = pd.read_csv(counts_path, sep=',')
    counts.columns = ['id', 'counts']
    counts = counts[[True if i in id_to_name else False for i in counts['id']]]
    counts = counts.assign(tf_name = [id_to_name[i] for i in counts['id']])
    counts = counts.set_index(['tf_name'])
    # counts


    # ANNOTATED ESDEG
    esdeg = esdeg.assign(counts = [counts.loc[i]['counts'] if i in counts.index else 0 for i in esdeg['tf_name']])
    esdeg = esdeg.assign(lfc = [deg.loc[i]['log2FoldChange'] if i in deg.index else 0 for i in esdeg['tf_name']])
    esdeg = esdeg.assign(de_padj = [deg.loc[i]['padj'] if i in deg.index else 0 for i in esdeg['tf_name']])

    if filter_flag:
        esdeg = esdeg[esdeg['me_padj'] < me_padj_thr]
        esdeg = esdeg[esdeg['log2(or)'] > lor_thr]
        summary_deg = esdeg[np.logical_and(esdeg['de_padj'] < de_padj_thr, np.abs(esdeg['lfc']) > lfc_thr)]
        summary_counts = esdeg[esdeg['counts'] > counts_filter]
        esdeg = pd.concat([summary_deg, summary_counts])
        esdeg = esdeg[~esdeg.index.duplicated(keep='first')]

    esdeg = esdeg.sort_values(by='de_padj')
    if best_flag:
        esdeg = esdeg.drop_duplicates(subset=['tf_name'], keep='first')

    if 'jaspar_cluster' in esdeg.columns:
        esdeg = esdeg[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'jaspar_cluster', 'log2(or)', 'me_padj', 'counts', 'lfc', 'de_padj', 'genes']]
    else:
        esdeg = esdeg[['motif_id', 'tf_name', 'tf_class', 'tf_family', 'log2(or)', 'me_padj', 'counts', 'lfc', 'de_padj', 'genes']]
        
    esdeg = esdeg.reset_index(drop=True)
    return esdeg
