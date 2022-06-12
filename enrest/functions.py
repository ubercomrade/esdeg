import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm, chi2
from itertools import tee
from enrest.speedup import montecarlo_enrichment, montecarlo_fraction

def get_threshold(scores):
    container = []
    number_of_sites = len(scores)
    last_score = scores[0]
    for count, score in enumerate(scores[1:], 1):
        if score == last_score:
            continue
        elif count/number_of_sites > 0.0005:
            container.append((last_score, count/number_of_sites))
            break
        elif score != last_score:
            container.append((last_score, count/number_of_sites))
            last_score = score 
    return container


def calculate_fprs(min_fpr, n=3):
    container = []
    step = (np.log10(np.min(min_fpr)) - np.log10(0.0005)) / n
    for i in range(n):
        fpr = 10**(np.log10(0.0005) + (i*step))
        container.append(fpr)
    container.append(min_fpr)
    return container


def split_scores_by_gene_ids(scores, gc, all_ids, deg_ids, other_ids):
    ids, index_all, index_deg = np.intersect1d(all_ids, deg_ids, assume_unique=True, return_indices=True) 
    deg_scores = scores[index_all]
    deg_gc = gc[index_all]
    genes = all_ids[index_all]
    
    ids, index_all, index_other = np.intersect1d(all_ids, other_ids, assume_unique=True, return_indices=True) 
    other_scores = scores[index_all]
    other_gc = gc[index_all]
    return deg_scores, deg_gc, other_scores, other_gc, genes


def hartung(p):
    '''
     Hartung, J. (1999): "A note on combining dependent tests of significance",
                         Biometrical Journal, 41(7), 849--855.
    '''
    L = np.ones(len(p), dtype=float) # zeros weight
    t = norm.ppf(p)
    n = float(len(p))
    avt = np.sum(t)/n
    q = np.sum((t - avt)**2)/(n-1)  # Hartung, eqn. (2.2)
    rhohat = 1 - q
    rhostar = max(-1/(n-1), rhohat) # Hartung, p. 851
    kappa = (1 + 1/(n-1) - rhostar)/10 # Hartung, p. 853
    Ht = np.sum(L*t)/np.sqrt(np.sum(L**2)+((np.sum(L))**2-np.sum(L**2))*(rhostar+kappa*np.sqrt(2/(n-1))*(1-rhostar))) # Hartung, p. 854, eq 2.4
    pvalue=norm.cdf(Ht)
    return pvalue


def get_deg_gene_ids(df, cond, padj_thr=0.05, log2fc_thr=1):
    log2fc_thr = abs(log2fc_thr)
    if cond == 'ALL':
        df = df[np.logical_or(df['log2FoldChange'] <= -log2fc_thr, df['log2FoldChange'] >= log2fc_thr)]
    elif cond == 'DOWN':
        df = df[df['log2FoldChange'] <= -log2fc_thr]
    elif cond == 'UP':
        df = df[df['log2FoldChange'] >= log2fc_thr]
    df = df[df['padj'] <= padj_thr]
    gene_ids = np.array([i for i in df['id'] if isinstance(i, str)])
    return gene_ids


def get_other_gene_ids_for_deg_case(df, padj_thr=0.05, log2fc_thr=np.log2(5/4)):
    log2fc_thr = abs(log2fc_thr)
    df = df[np.logical_and(df['log2FoldChange'] >= -log2fc_thr, df['log2FoldChange'] <= log2fc_thr)]
    df = df[df['padj'] > padj_thr]
    gene_ids = np.array([i for i in df['id'] if isinstance(i, str)])
    return gene_ids


def get_other_gene_ids_for_set_case(set_ids, all_ids):
    gene_ids = set(list(all_ids)) - set(list(set_ids))
    return np.array(list(gene_ids))


def run_test(genes, set_scores, set_gc, other_scores, other_gc, threshold_table, gc_threshold, parameter):
    results = {('LOW', 'log2Fold'): 0, ('LOW', 'pval'): 0, ('LOW', 'genes'): [],
               ('MIDDLE', 'log2Fold'): 0, ('MIDDLE', 'pval'): 0, ('MIDDLE', 'genes'): [],
               ('HIGH', 'log2Fold'): 0, ('HIGH', 'pval'): 0, ('HIGH', 'genes'): [],
               ('COMBINE', 'pval'): 0}
    pvals = []
    genes = np.array(genes)
    for index, level in zip(range(0, len(threshold_table)), ['LOW', 'MIDDLE', 'HIGH']):
        threshold, fpr = threshold_table[index]
        if parameter == "enrichment":
            z_score, fold = montecarlo_enrichment(set_scores, set_gc,
                                                  other_scores, other_gc,
                                                  threshold, gc_threshold)
            genes_with_bs = genes[np.sum(np.greater_equal(set_scores, threshold), axis=1) > 0]
        elif parameter == "fraction":
            z_score, fold = montecarlo_fraction(set_scores, set_gc,
                                                other_scores, other_gc,
                                                threshold, gc_threshold)
            genes_with_bs = genes[np.greater_equal(set_scores, threshold)]
        pval = stats.norm.sf(z_score)          
        pvals.append(pval)
        results[(level, 'pval')] = pval
        results[(level, 'log2Fold')] = fold
        results[(level, 'genes')] = ';'.join(genes_with_bs)
    cpval = hartung(np.array(pvals))
    results[('COMBINE', 'pval')] = cpval
    return results


def run_test_fasta(set_scores, set_gc, other_scores, other_gc, threshold_table, gc_threshold, parameter):
    results = {('LOW', 'log2Fold'): 0, ('LOW', 'pval'): 0,
               ('MIDDLE', 'log2Fold'): 0, ('MIDDLE', 'pval'): 0,
               ('HIGH', 'log2Fold'): 0, ('HIGH', 'pval'): 0,
               ('COMBINE', 'pval'): 0}
    pvals = []
    for index, level in zip(range(0, len(threshold_table)), ['LOW', 'MIDDLE', 'HIGH']):
        threshold, fpr = threshold_table[index]
        if parameter == "enrichment":
            z_score, fold = montecarlo_enrichment(set_scores, set_gc,
                                                  other_scores, other_gc,
                                                  threshold, gc_threshold)
        elif parameter == "fraction":
            z_score, fold = montecarlo_fraction(set_scores, set_gc,
                                                other_scores, other_gc,
                                                threshold, gc_threshold)
        pval = stats.norm.sf(z_score)          
        pvals.append(pval)
        results[(level, 'pval')] = pval
        results[(level, 'log2Fold')] = fold
    cpval = hartung(np.array(pvals))
    results[('COMBINE', 'pval')] = cpval
    return results


def write_table(head, data, path):
    with open(path, 'w') as file:
        file.write(head)
        for line in data:
            file.write('\t'.join(map(str, line)) + '\n')
    pass