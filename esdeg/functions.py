import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm, chi2
from itertools import tee


def calculate_gc(sequences):
    length = sequences.shape[1]
    gc = np.sum(np.logical_or(sequences == 1, sequences == 2), axis=1)
    gc = gc / length
    return gc


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


def run_test(genes, foreground, foreground_gc, other, other_gc, gc_threshold, parameter):    
    if parameter == "fraction":
        #from number of sites to number of peaks with site
        foreground = np.array(foreground >= 1, dtype=int)
        other = np.array(other >= 1, dtype=int)
        
    #run montecarlo
    if gc_threshold < 1.0:
        z_scores, odds_ratio = montecarlo_gc(foreground, foreground_gc, 
                                          other, other_gc, gc_threshold)
    else:
        z_scores, odds_ratio = montecarlo(foreground, other)
    
    #list of genes with site (thr 0.0005)
    genes_with_bs = genes[np.greater_equal(foreground[:,1], 1)]
    
    #calculate one -log10(p-value) from several p-values by using Hartung method
    pv = hartung(stats.norm.sf(z_scores))
    lpv = -np.log10(pv)
    
    #calculate distance
    if np.log2(odds_ratio) > 0:
        distance = np.sqrt(np.power(odds_ratio, 2) + np.power(lpv, 2))
    else:
        distance = -np.sqrt(np.power(1 / odds_ratio, 2) + np.power(lpv, 2))
        
    #write results
    results = {'log(or)': np.log2(odds_ratio),
               'distance': distance,
               'pval': pv,
               'genes': ';'.join(list(genes_with_bs))}
    return results


def montecarlo_gc(foreground, foreground_gc, other, other_gc, gc_threshold):
    number_of_foreground, number_of_thresholds = foreground.shape
    number_of_other = other.shape[0]
    gc_index = np.abs([other_gc - i for i in foreground_gc])
    gc_index = gc_index < gc_threshold
    gc_index = [np.where(i)[0] for i in gc_index]
    random = np.zeros((1000, number_of_thresholds), dtype=np.int64)
    real = np.sum(foreground, axis=0)
    for i in range(1000):
        index = [np.random.choice(j) for j in gc_index]
        sample = other[index]
        random[i] = np.sum(sample, axis=0)
    random_mean, random_std = np.mean(random, axis=0), np.std(random, axis=0)
    z_score = np.abs((real - random_mean) / random_std)
    return z_score, np.mean(real / np.mean(random, axis=0))


def montecarlo(foreground, other):
    number_of_foreground, number_of_thresholds = foreground.shape
    number_of_other = other.shape[0]
    random = np.zeros((1000, number_of_thresholds), dtype=np.int64)
    real = np.sum(foreground, axis=0)
    for i in range(1000):
        index = np.random.choice(number_of_other, number_of_foreground)
        sample = other[index]
        random[i] = np.sum(sample, axis=0)
    random_mean, random_std = np.mean(random, axis=0), np.std(random, axis=0)
    z_score = np.abs((real - random_mean) / random_std)
    return z_score, np.mean(real / np.mean(random, axis=0))


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


def split_by_gene_ids(counts, gc, ids, foreground_ids, other_ids):
    _, index, index_foreground = np.intersect1d(ids, foreground_ids, assume_unique=True, return_indices=True) 
    foreground_counts = counts[index]
    foreground_gc = gc[index]
    genes = ids[index]
    
    _, index, index_other = np.intersect1d(ids, other_ids, assume_unique=True, return_indices=True) 
    other_counts = counts[index]
    other_gc = gc[index]
    return foreground_counts, foreground_gc, other_counts, other_gc, genes


def get_deg_gene_ids(df, cond, padj_thr=0.05, log2fc_thr=1):
    log2fc_thr = abs(log2fc_thr)
    if cond == 'all':
        df = df[np.logical_or(df['log2FoldChange'] <= -log2fc_thr, df['log2FoldChange'] >= log2fc_thr)]
    elif cond == 'down':
        df = df[df['log2FoldChange'] <= -log2fc_thr]
    elif cond == 'up':
        df = df[df['log2FoldChange'] >= log2fc_thr]
    df = df[df['padj'] <= padj_thr]
    gene_ids = np.array([i for i in df['id'] if isinstance(i, str)])
    return np.unique(gene_ids)


def get_other_gene_ids_for_deg_case(df, padj_thr=0.05, log2fc_thr=np.log2(5/4)):
    log2fc_thr = abs(log2fc_thr)
    #df = df[np.logical_and(df['log2FoldChange'] >= -log2fc_thr, df['log2FoldChange'] <= log2fc_thr)]
    df = df[df['padj'] > padj_thr]
    gene_ids = np.array([i for i in df['id'] if isinstance(i, str)])
    return np.unique(gene_ids)


def get_other_gene_ids_for_set_case(set_ids, all_ids):
    gene_ids = set(list(all_ids)) - set(list(set_ids))
    return np.array(list(gene_ids))