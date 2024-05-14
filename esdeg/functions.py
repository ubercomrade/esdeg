import os
import json
import numpy as np
import pandas as pd
import collections
from itertools import product
from scipy import stats
from scipy.stats import norm, chi2
from itertools import tee
from numpy.lib.stride_tricks import sliding_window_view


#str to int
def seq_to_int(sequences):
    number_of_promoters = len(sequences)
    length_of_promoters = len(sequences[0])
    converter = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    results = np.zeros((number_of_promoters, length_of_promoters), dtype=int)
    for index in range(number_of_promoters):
        results[index] = np.array([converter[i] for i in sequences[index]])
    return(results)

#PWM
def dict_to_array(motif):
    motif = [motif[i] for i in motif.keys()]
    return np.array(motif)


def pfm_to_pwm(pfm):
    background = 0.25
    pwm = np.log2(pfm / background)
    return pwm


def pfm_to_pwm_cisbp(pfm):
    background = 0.25
    pwm = np.log2((pfm + 0.01) / background)
    return pwm


def pcm_to_pfm(pcm):
    number_of_sites = pcm.sum(axis=0)
    nuc_pseudo = 0.25
    pfm = (pcm + nuc_pseudo) / (number_of_sites + 1)
    return pfm


def min_score(matrix):
    return np.sum(matrix.min(axis=0))


def max_score(matrix):
    return np.sum(matrix.max(axis=0))


def to_score(matrix, norm_value):
    min_s = min_score(matrix)
    max_s = max_score(matrix)
    score = norm_value * (max_s - min_s) + min_s
    return score


def pfm_to_pwm_hocomoco(pfm):
    background = 0.25
    pwm = np.log2((pfm + 0.01) / background)
    return pwm


def pfm_to_pwm_cisbp(pfm):
    background = 0.25
    pwm = np.log2((pfm + 0.01) / background)
    return pwm


def hocomoco_to_pwm(motif_data):
    pfm = np.array(motif_data['pfm']).T
    pwm = pfm_to_pwm_hocomoco(pfm)
    return pwm


def jaspar_to_pwm(motif):
    pcm = motif.counts
    pcm = dict_to_array(pcm)
    pfm = pcm_to_pfm(pcm)
    pwm = pfm_to_pwm(pfm)
    return pwm


def scores_of_sites(pwm, matrix_of_sites):
    return np.sum(np.take_along_axis(pwm, matrix_of_sites, axis=0), axis=1)


def all_scores_pwm(numeric_sequenes, pwm):
    length_of_site = pwm.shape[1]
    number_of_promoters = numeric_sequenes.shape[0]
    number_of_scores = numeric_sequenes.shape[1] - length_of_site + 1
    scores = np.zeros((number_of_promoters, number_of_scores * 2), dtype=float)
    for i in range(number_of_promoters):
        matrix_of_sites = sliding_window_view(numeric_sequenes[i], length_of_site)
        matrix_of_sites_rc = np.flip(np.array((matrix_of_sites -3) * -1), axis=1)  # reverse+copliment
        scores_forward = scores_of_sites(pwm, matrix_of_sites)
        scores_backward = scores_of_sites(pwm, matrix_of_sites_rc)
        scores[i] = np.concatenate((scores_forward, scores_backward))
    return scores


def count_sites(scores, threshold_table):
    counts = np.zeros((scores.shape[0], threshold_table.shape[0]), dtype=int)
    for i in range(threshold_table.shape[0]):
        counts[:,i] = np.sum(np.greater_equal(scores, threshold_table[i][0]), axis=1)
    return counts

# def calculate_gc(sequences):
#     length = sequences.shape[1]
#     gc = np.sum(np.logical_or(sequences == 1, sequences == 2), axis=1)
#     gc = gc / length
#     return gc


def calculate_gc(sequences):
    gc = []
    length = len(sequences[0])
    for seq in sequences:
        counter = collections.Counter(seq)
        gc.append((counter['C'] + counter['G']) / length)
    return np.array(gc)


# def make_k_mers(order):
#     tmp = product('ACGT', repeat=order + 1)
#     k_mer = []
#     for i in tmp:
#         k_mer.append(''.join(i))
#     k_mer_dict = dict()
#     index = 0
#     for i in k_mer:
#         k_mer_dict[i] = index
#         index += 1
#     return(k_mer_dict)


def get_threshold(scores):
    container = []
    number_of_sites = len(scores)
    last_score = scores[0]
    for count, score in enumerate(scores[1:], 1):
        if score == last_score:
            continue
        elif count/number_of_sites > 0.001:
            container.append((last_score, count/number_of_sites))
            break
        elif score != last_score:
            container.append((last_score, count/number_of_sites))
            last_score = score
    return container


# def calculate_fprs(min_fpr, n=3):
#     container = []
#     step = (np.log10(np.min(min_fpr)) - np.log10(0.0005)) / n
#     for i in range(n):
#         fpr = 10**(np.log10(0.0005) + (i*step))
#         container.append(fpr)
#     container.append(min_fpr)
#     return container


def run_test(genes, foreground, foreground_gc, other, other_gc, gc_threshold, parameter, number_of_uniq_fprs):
    if parameter == "fraction":
        #from number of sites to number of peaks with site
        foreground = np.array(foreground >= 1, dtype=int)
        other = np.array(other >= 1, dtype=int)

    #run montecarlo
    z_scores, odds_ratios = montecarlo(foreground, foreground_gc,
                                       other, other_gc,
                                       gc_threshold)
    z_scores = -1.0 * z_scores # for `one sided` and `greater`
    log_pvalues=norm.logcdf(z_scores) # natural log
    log_pvalue = bonferroni_combined_log_pvalues(log_pvalues, number_of_uniq_fprs)
    #log_pvalue = hartung_log(z_scores)
    index_of_best = np.argmin(log_pvalues)

    genes_best = genes[np.greater_equal(foreground[:,index_of_best], 1)]
    odds_ratio_best = odds_ratios[index_of_best]

    #write results
    results = {'log2(or)': np.log2(np.max(odds_ratios)),
               'ln(pval)': log_pvalue,
               'genes': ';'.join(list(genes_best))}
    return results, log_pvalues, odds_ratios, index_of_best


# def run_test(genes, foreground, foreground_gc, other, other_gc, gc_threshold, parameter):
#     if parameter == "fraction":
#         #from number of sites to number of peaks with site
#         foreground = np.array(foreground >= 1, dtype=int)
#         other = np.array(other >= 1, dtype=int)
#
#     #run montecarlo
#     z_scores, odds_ratio = montecarlo(foreground, foreground_gc,
#                                           other, other_gc, gc_threshold)
#     #z_scores, odds_ratio = montecarlo(foreground, other)
#
#     #list of genes with site (thr 0.0005)
#     genes_low_thr = genes[np.greater_equal(foreground[:,0], 1)]
#
#     #list of genes with site (thr 0.0001)
#     genes_high_thr = genes[np.greater_equal(foreground[:,-1], 1)]
#
#
#     #calculate one ln(p-value) from several p-values by using Hartung method
#     #pv = hartung(stats.norm.sf(z_scores)) # old. work with p-values
#     lpv = hartung_log(z_scores) # new. work with ln(p-values)
#     #lpv = -np.log10(pv)
#
# #    #calculate distance
# #     if np.log2(odds_ratio) > 0:
# #         distance = np.sqrt(np.power(odds_ratio, 2) + np.power(lpv, 2))
# #     else:
# #         distance = -np.sqrt(np.power(1 / odds_ratio, 2) + np.power(lpv, 2))
#
#     #write results
#     results = {'log2(or)': np.log2(odds_ratio),
#                'ln(pval)': lpv,
#                'genes_low_thr': ';'.join(list(genes_low_thr)),
#                'genes_high_thr': ';'.join(list(genes_high_thr))}
#     return results


def montecarlo(foreground, foreground_gc, other, other_gc, gc_threshold):
    number_of_foreground, number_of_thresholds = foreground.shape
    number_of_other = other.shape[0]
    gc_index = np.abs([other_gc - i for i in foreground_gc])
    gc_index = gc_index < gc_threshold
    number_of_iterations = 1000
    indexes = np.zeros((number_of_foreground, number_of_iterations), dtype=np.int64)
    for i in range(number_of_foreground):
        indexes[i] = np.random.choice(gc_index[i].nonzero()[0], number_of_iterations)
    random = np.zeros((number_of_iterations, number_of_thresholds), dtype=np.int64)
    real = np.sum(foreground, axis=0)
    for i in range(number_of_iterations):
        index = indexes[:,i]
        sample = other[index]
        random[i] = np.sum(sample, axis=0)
    random_mean, random_std = np.mean(random, axis=0), np.std(random, axis=0)
    z_score = (real - random_mean) / random_std
    return z_score, real / random_mean


# def montecarlo(foreground, other):
#     number_of_foreground, number_of_thresholds = foreground.shape
#     number_of_other = other.shape[0]
#     random = np.zeros((1000, number_of_thresholds), dtype=np.int64)
#     real = np.sum(foreground, axis=0)
#     for i in range(1000):
#         index = np.random.choice(number_of_other, number_of_foreground)
#         sample = other[index]
#         random[i] = np.sum(sample, axis=0)
#     random_mean, random_std = np.mean(random, axis=0), np.std(random, axis=0)
#     z_score = np.abs((real - random_mean) / random_std)
#     return z_score, np.mean(real / np.mean(random, axis=0))


#COMBINED PVALUES
def bonferroni_combined_pvalues(pvals, number_of_uniq_fprs):
    return np.min([1.0, np.min(pvals) * number_of_uniq_fprs])

def bonferroni_combined_log_pvalues(log_pvals, number_of_uniq_fprs):
    return np.min([0.0, np.min(log_pvals) + np.log(number_of_uniq_fprs)])

# def hartung(p):
#     '''
#      Hartung, J. (1999): "A note on combining dependent tests of significance",
#                          Biometrical Journal, 41(7), 849--855.
#     '''
#     L = np.ones(len(p), dtype=float) # zeros weight
#     t = norm.ppf(p)
#     n = float(len(p))
#     avt = np.sum(t)/n
#     q = np.sum((t - avt)**2)/(n-1)  # Hartung, eqn. (2.2)
#     rhohat = 1 - q
#     rhostar = max(-1/(n-1), rhohat) # Hartung, p. 851
#     kappa = (1 + 1/(n-1) - rhostar)/10 # Hartung, p. 853
#     Ht = np.sum(L*t)/np.sqrt(np.sum(L**2)+((np.sum(L))**2-np.sum(L**2))*(rhostar+kappa*np.sqrt(2/(n-1))*(1-rhostar))) # Hartung, p. 854, eq 2.4
#     pvalue=norm.cdf(Ht)
#     return pvalue

def hartung_log(z_scores):
    '''
     Hartung, J. (1999): "A note on combining dependent tests of significance",
                         Biometrical Journal, 41(7), 849--855.
    '''
    L = np.ones(len(z_scores)) # zeros weight
    n = float(len(z_scores))
    avt = np.sum(z_scores)/n
    q = np.sum((z_scores - avt)**2)/(n-1)  # Hartung, eqn. (2.2)
    rhohat = 1 - q
    rhostar = max(-1/(n-1), rhohat) # Hartung, p. 851
    #kappa = (1 + 1/(n-1) - rhostar)/10 # Hartung, p. 853
    kappa = 0.2
    combined_z_score = np.sum(L*z_scores)/np.sqrt(np.sum(L**2)+((np.sum(L))**2-np.sum(L**2))*(rhostar+kappa*np.sqrt(2/(n+1))*(1-rhostar))) # Hartung, p. 854, eq 2.4
    log_pvalue = norm.logcdf(combined_z_score)
    return log_pvalue

#calculate adj.pvalues.
# def fdr(p_vals):
#     ranked_p_values = stats.rankdata(p_vals)
#     fdr = p_vals * len(p_vals) / ranked_p_values
#     fdr[fdr > 1] = 1
#     return fdr
#
# def fdr_log(lg_pvals):
#     ranked_p_values = stats.rankdata(lg_pvals)
#     fdr = lg_pvals + np.log(len(lg_pvals)) - np.log(ranked_p_values)
#     fdr[fdr > 0] = 0
#     return fdr


#calculate adj.pvalues. Based on function statsmodels.stats.multitest.fdrcorrection from statsmodels
def _ecdf_log(x):
    nobs = len(x)
    return np.log(np.arange(1,nobs+1)) - np.log(nobs)

def fdrcorrection_log(lg_pvals):
    pvals = np.asarray(lg_pvals)
    lg_pvals_sortind = np.argsort(lg_pvals)
    lg_pvals_sorted = np.take(pvals, lg_pvals_sortind)

    ecdffactor = _ecdf_log(lg_pvals_sorted)

    lg_pvals_corrected_raw = lg_pvals_sorted - ecdffactor
    lg_pvals_corrected = np.minimum.accumulate(lg_pvals_corrected_raw[::-1])[::-1]
    del lg_pvals_corrected_raw
    lg_pvals_corrected[lg_pvals_corrected > 0] = 0
    lg_pvals_corrected_ = np.empty_like(lg_pvals_corrected)
    lg_pvals_corrected_[lg_pvals_sortind] = lg_pvals_corrected
    del lg_pvals_corrected
    return lg_pvals_corrected_


#work_with genes
def split_by_gene_ids(counts, gc, ids, foreground_ids, other_ids):
    _, index, index_foreground = np.intersect1d(ids, foreground_ids, assume_unique=False, return_indices=True)
    foreground_counts = counts[index]
    foreground_gc = gc[index]
    genes = ids[index]

    _, index, index_other = np.intersect1d(ids, other_ids, assume_unique=False, return_indices=True)
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
    df = df[np.logical_and(df['log2FoldChange'] >= -log2fc_thr, df['log2FoldChange'] <= log2fc_thr)]
    df = df[df['padj'] > padj_thr]
    gene_ids = np.array([i for i in df['id'] if isinstance(i, str)])
    return np.unique(gene_ids)


def get_other_gene_ids_for_set_case(set_ids, all_ids):
    gene_ids = set(list(all_ids)) - set(list(set_ids))
    return np.array(list(gene_ids))


#clusters
def get_motif_to_cluster(cluster_path):
    clusters = pd.read_csv(cluster_path, sep='\t')
    motif_to_cluster = dict()
    for index, line in clusters.iterrows():
        cluster_id = line[0]
        for motif_id in line['motif_ids'].split(','):
            motif_to_cluster[motif_id] = cluster_id
    return motif_to_cluster
