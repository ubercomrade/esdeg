#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy import stats
from numba import njit, float64, int64


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


def choose_threshold(path, fpr_for_thr):
    container = list()
    append = container.append
    with open(path, 'r') as file:
        for line in file:
            append(tuple(map(float, line.strip().split())))
    file.close()
    container = sorted(container, key=itemgetter(1))
    last_score, last_fpr = container[0]
    for line in container:
        if line[1] > fpr_for_thr:
            break
        else:
            last_score, last_fpr = line
    return(last_score)


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


@njit
def calculate_enrichment(scores, threshold):
    counts = np.sum(np.greater_equal(scores, threshold))
    total = scores.shape[0] * scores.shape[1]
    enrichment = counts / total
    if enrichment == 0:
        enrichment = 10**(-6)
    return enrichment


@njit
def montecarlo(deg_scores, other_scores, threshold):
    real_enrichmnet = calculate_enrichment(deg_scores, threshold)
    number_of_deg = len(deg_scores)
    indexes = np.arange(len(other_scores))
    vec_random_enrichment = np.zeros(1000)
    for i in range(1000):
        sample_indexes = np.random.choice(indexes, number_of_deg)
        sample = other_scores[sample_indexes]
        vec_random_enrichment[i] = calculate_enrichment(sample, threshold)
    random_enrichmnet_mean, random_enrichmnet_std = np.mean(vec_random_enrichment), np.std(vec_random_enrichment)
    z_score = (real_enrichmnet - random_enrichmnet_mean) / random_enrichmnet_std
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