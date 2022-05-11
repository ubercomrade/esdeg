import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm, chi2
from numba import njit, float64, int64
from itertools import tee


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


def split_scores_by_gene_ids(scores, all_ids, deg_ids, other_ids):
    deg_scores = []
    other_scores = []
    exist = []
    set_deg_ids = set(deg_ids)
    set_other_ids = set(other_ids)
    for index, i in enumerate(all_ids):
        if i in set_deg_ids:
            exist.append(i)
            deg_scores.append(scores[index])
        elif i in set_other_ids:
            other_scores.append(scores[index])
        else:
            continue
    return np.array(deg_scores), np.array(other_scores), exist


@njit(cache=True)
def calculate_enrichment(scores, threshold):
    counts = np.sum(np.greater_equal(scores, threshold))
    total = scores.shape[0] * scores.shape[1]
    enrichment = counts / total
    if enrichment == 0.0:
        enrichment = 10**(-6)
    return enrichment


@njit(cache=True)
def montecarlo_enrichment(deg_scores, other_scores, threshold):
    real_enrichmnet = calculate_enrichment(deg_scores, threshold)
    number_of_deg = len(deg_scores)
    indexes = np.arange(len(other_scores))
    vec_random_enrichment = np.zeros(1000)
    for i in range(1000):
        sample_indexes = np.random.choice(indexes, number_of_deg)
        sample = other_scores[sample_indexes]
        vec_random_enrichment[i] = calculate_enrichment(sample, threshold)
    random_enrichmnet_mean, random_enrichmnet_std = np.mean(vec_random_enrichment), np.std(vec_random_enrichment)
    z_score = abs((real_enrichmnet - random_enrichmnet_mean) / random_enrichmnet_std)
    #print(real_enrichmnet, random_enrichmnet_mean, real_enrichmnet / random_enrichmnet_mean)
    return z_score, np.log2(real_enrichmnet / random_enrichmnet_mean)


@njit(cache=True)
def montecarlo_fraction(deg_scores, other_scores, threshold):
    number_of_deg = len(deg_scores)
    number_of_deg_with_tfbs = np.sum(np.greater_equal(deg_scores, threshold))
    real_fraction = number_of_deg_with_tfbs / number_of_deg
    indexes = np.arange(len(other_scores))
    vec_random_fraction = np.zeros(1000, dtype=np.float64)
    for i in range(1000):
        sample_indexes = np.random.choice(indexes, number_of_deg)
        sample = other_scores[sample_indexes]
        number_of_other_with_tfbs = np.sum(np.greater_equal(sample, threshold))
        vec_random_fraction[i] = number_of_other_with_tfbs / number_of_deg
    random_fraction_mean, random_fraction_std = np.mean(vec_random_fraction), np.std(vec_random_fraction)
    z_score = abs((real_fraction - random_fraction_mean) / random_fraction_std)
    #print(real_fraction, random_fraction_mean, real_fraction / random_fraction_mean)
    return z_score, np.log2(real_fraction / random_fraction_mean)


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
    gene_ids = [i for i in df['id'] if isinstance(i, str)]
    return gene_ids


def get_other_gene_ids_for_deg_case(df, padj_thr=0.05, log2fc_thr=np.log2(5/4)):
    log2fc_thr = abs(log2fc_thr)
    #df = df[np.logical_and(df['log2FoldChange'] >= -log2fc_thr, df['log2FoldChange'] <= log2fc_thr)]
    df = df[df['padj'] > padj_thr]
    gene_ids = [i for i in df['id'] if isinstance(i, str)]
    return gene_ids


def get_other_gene_ids_for_set_case(set_ids, all_ids):
    gene_ids = all_ids - set_ids
    return gene_ids


def run_test(genes, set_scores, other_scores, threshold_table, parameter):
    results = {('LOW', 'log2Fold'): 0, ('LOW', 'pval'): 0, ('LOW', 'genes'): [],
               ('MIDDLE', 'log2Fold'): 0, ('MIDDLE', 'pval'): 0, ('MIDDLE', 'genes'): [],
               ('HIGH', 'log2Fold'): 0, ('HIGH', 'pval'): 0, ('HIGH', 'genes'): [],
               ('COMBINE', 'pval'): 0}
    pvals = []
    genes = np.array(genes)
    for index, level in zip(range(0, len(threshold_table)), ['LOW', 'MIDDLE', 'HIGH']):
        threshold, fpr = threshold_table[index]
        if parameter == "enrichment":
            z_score, fold = montecarlo_enrichment(set_scores, other_scores, threshold)
            genes_with_bs = genes[np.sum(np.greater_equal(set_scores, threshold), axis=1) > 0]
        elif parameter == "fraction":
            z_score, fold = montecarlo_fraction(set_scores, other_scores, threshold)
            genes_with_bs = genes[np.greater_equal(set_scores, threshold)]
        pval = stats.norm.sf(z_score)          
        pvals.append(pval)
        results[(level, 'pval')] = pval
        results[(level, 'log2Fold')] = fold
        results[(level, 'genes')] = ';'.join(genes_with_bs)
    cpval = hartung(np.array(pvals))
    results[('COMBINE', 'pval')] = cpval
    return results


def run_test_fasta(set_scores, other_scores, threshold_table, parameter):
    results = {('LOW', 'log2Fold'): 0, ('LOW', 'pval'): 0,
               ('MIDDLE', 'log2Fold'): 0, ('MIDDLE', 'pval'): 0,
               ('HIGH', 'log2Fold'): 0, ('HIGH', 'pval'): 0,
               ('COMBINE', 'pval'): 0}
    pvals = []
    for index, level in zip(range(0, len(threshold_table)), ['LOW', 'MIDDLE', 'HIGH']):
        threshold, fpr = threshold_table[index]
        if parameter == "enrichment":
            z_score, fold = montecarlo_enrichment(set_scores, other_scores, threshold)
        elif parameter == "fraction":
            z_score, fold = montecarlo_fraction(set_scores, other_scores, threshold)
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
