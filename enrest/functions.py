#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats
from numba import njit, float64, int64
from itertools import tee


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)
    

def get_best_scores(scan):
    container = []
    for gname, s in scan:
        container.append((gname, s.max()))
    return container


def get_all_scores(scan):
    all_scores = np.zeros((len(scan), len(scan[0][1])))
    for index, line in enumerate(scan):
        all_scores[index] = line[1]
    return all_scores


def get_all_flatten_scores(scan):
    all_scores = get_all_scores(scan)
    all_scores = all_scores.ravel()
    all_scores.sort()
    all_scores = all_scores[::-1]
    return all_scores


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


# def calculate_fpr(n):
#     if n == 0:
#         return(0.0005)
#     else:
#         parameter = 1.5
#         fpr = 0.0005
#         for i in range(n):
#             fpr = fpr / parameter
#             parameter = 1 + parameter / 4 # old val = 2
#         return(fpr)


def calculate_fprs(min_fpr, n=3):
    container = []
    step = (np.log10(np.min(min_fpr)) - np.log10(0.0005)) / n
    for i in range(n):
        fpr = 10**(np.log10(0.0005) + (i*step))
        container.append(fpr)
    container.append(min_fpr)
    return container


def split_scores_by_gene_ids(scan_results, deg_ids, other_ids):
    deg_scores = []
    other_scores = []
    exist = []
    for gname, score in scan_results:
        if gname in deg_ids:
            exist.append(gname)
            deg_scores.append(score)
        elif gname in other_ids:
            other_scores.append(score)
        else:
            continue
    #print(set(deg_ids) - set(exist))
    return np.array(deg_scores), np.array(other_scores)


def fisher_test(deg_scores, other_scores, threshold):
    number_of_deg_with_tfbs = np.sum(np.greater_equal(deg_scores, threshold))
    number_of_deg_without_tfbs = len(deg_scores) - number_of_deg_with_tfbs
    number_of_other_with_tfbs = np.sum(np.greater_equal(other_scores, threshold))
    number_of_other_without_tfbs = len(other_scores) - number_of_other_with_tfbs
    table = [[number_of_deg_with_tfbs, number_of_other_with_tfbs],
        [number_of_deg_without_tfbs, number_of_other_without_tfbs]]
    oddsr, fisher_pval = stats.fisher_exact(table, alternative='two-sided') # {‘two-sided’, ‘less’, ‘greater’}
    return fisher_pval


# @njit(cache=True)
# def calculate_enrichment(scores, threshold):
#     counts = np.sum(np.greater_equal(scores, threshold))
#     total = scores.shape[0] * scores.shape[1]
#     enrichment = counts / total
#     if enrichment == 0:
#         enrichment = 10**(-6)
#     return enrichment


# @njit(cache=True)
# def montecarlo_enrichment(deg_scores, other_scores, threshold):
#     real_enrichmnet = calculate_enrichment(deg_scores, threshold)
#     number_of_deg = len(deg_scores)
#     indexes = np.arange(len(other_scores))
#     vec_random_enrichment = np.zeros(1000)
#     for i in range(1000):
#         sample_indexes = np.random.choice(indexes, number_of_deg)
#         sample = other_scores[sample_indexes]
#         vec_random_enrichment[i] = calculate_enrichment(sample, threshold)
#     random_enrichmnet_mean, random_enrichmnet_std = np.mean(vec_random_enrichment), np.std(vec_random_enrichment)
#     z_score = (real_enrichmnet - random_enrichmnet_mean) / random_enrichmnet_std
#     return z_score


# @njit(cache=True)
# def montecarlo_fraction(deg_scores, other_scores, threshold):
#     number_of_deg = len(deg_scores)
#     number_of_deg_with_tfbs = np.sum(np.greater_equal(deg_scores, threshold))
#     real_fraction = number_of_deg_with_tfbs / number_of_deg
#     indexes = np.arange(len(other_scores))
#     vec_random_fraction = np.zeros(1000, dtype=np.float64)
#     for i in range(1000):
#         sample_indexes = np.random.choice(indexes, number_of_deg)
#         sample = other_scores[sample_indexes]
#         number_of_other_with_tfbs = np.sum(np.greater_equal(sample, threshold))
#         vec_random_fraction[i] = number_of_other_with_tfbs / number_of_deg
#     random_fraction_mean, random_fraction_std = np.mean(vec_random_fraction), np.std(vec_random_fraction)
#     z_score = (real_fraction - random_fraction_mean) / random_fraction_std
#     return z_score



@njit(fastmath=True)
def calculate_enrichment(scores, threshold_min, threshold_max):
    counts = np.sum(np.logical_and(np.greater_equal(scores, threshold_min),
                                   np.less(scores, threshold_max)))
    total = scores.shape[0] * scores.shape[1]
    enrichment = counts / total
    if enrichment == 0:
        enrichment = 10**(-6)
    return enrichment


@njit(cache=True, fastmath=False, parallel=False)
def montecarlo_enrichment(deg_scores, other_scores, threshold_min, threshold_max):
    real_enrichmnet = calculate_enrichment(deg_scores, threshold_min, threshold_max)
    number_of_deg = len(deg_scores)
    indexes = np.arange(len(other_scores))
    vec_random_enrichment = np.zeros(1000)
    for i in prange(1000):
        sample_indexes = np.random.choice(indexes, number_of_deg)
        sample = other_scores[sample_indexes]
        vec_random_enrichment[i] = calculate_enrichment(sample, threshold_min, threshold_max)
    random_enrichmnet_mean, random_enrichmnet_std = np.mean(vec_random_enrichment), np.std(vec_random_enrichment)
    z_score = (real_enrichmnet - random_enrichmnet_mean) / random_enrichmnet_std
    return z_score

    
@njit(cache=True, fastmath=False, parallel=False)
def montecarlo_fraction(deg_scores, other_scores, threshold_min, threshold_max):
    number_of_deg = len(deg_scores)
    number_of_deg_with_tfbs = np.sum(np.logical_and(np.greater_equal(deg_scores, threshold_min),
                                                   np.less(deg_scores, threshold_max)))
    real_fraction = number_of_deg_with_tfbs / number_of_deg
    indexes = np.arange(len(other_scores))
    vec_random_fraction = np.zeros(1000, dtype=np.float64)
    for i in prange(1000):
        sample_indexes = np.random.choice(indexes, number_of_deg)
        sample = other_scores[sample_indexes]
        number_of_other_with_tfbs = np.sum(np.logical_and(np.greater_equal(sample, threshold_min),
                                                   np.less(sample, threshold_max)))
        vec_random_fraction[i] = number_of_other_with_tfbs / number_of_deg
    random_fraction_mean, random_fraction_std = np.mean(vec_random_fraction), np.std(vec_random_fraction)
    z_score = (real_fraction - random_fraction_mean) / random_fraction_std
    return z_score


def get_deg_gene_ids(df, cond, padj_thr=0.1):
    if cond == 'ALL':
        df = df[np.logical_or(df['log2FoldChange'] <= -1, df['log2FoldChange'] >= 1)]
    elif cond == 'DOWN':
        df = df[df['log2FoldChange'] <= -1]
    elif cond == 'UP':
        df = df[df['log2FoldChange'] >= 1]
    df = df[df['padj'] <= padj_thr]
    gene_ids = [i for i in df['geneName'] if isinstance(i, str)]
    return gene_ids


def get_other_gene_ids(df, padj_thr=0.1):
    df = df[np.logical_and(df['log2FoldChange'] >= np.log2(4/5), df['log2FoldChange'] <= np.log2(5/4))]
    df = df[df['padj'] > padj_thr]
    gene_ids = [i for i in df['geneName'] if isinstance(i, str)]
    return gene_ids


def write_table(head, data, path):
    with open(path, 'w') as file:
        file.write(head)
        for line in data:
            file.write('\t'.join(map(str, line)) + '\n')
    pass