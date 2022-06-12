#%%pythran -Ofast -march=native -DUSE_XSIMD
import numpy as np
from numpy import random


#pythran export sort(float[])
def sort(scores):
    return np.sort(scores)[::-1]
    

#pythran export seq_to_int(str list)
def seq_to_int(sequences):
    number_of_promoters = len(sequences)
    length_of_promoters = len(sequences[0])
    vec = np.zeros((number_of_promoters, length_of_promoters), dtype=np.int64)
    for i in range(number_of_promoters):
        for j in range(length_of_promoters):
            if sequences[i][j] == 'A':
                vec[i][j] = 0
            elif sequences[i][j] == 'C':
                vec[i][j] = 1
            elif sequences[i][j] == 'G':
                vec[i][j] = 2
            elif sequences[i][j] == 'T':
                vec[i][j] = 3
            else:
                vec[i][j] = 4
    return vec
    

#pythran export get_score_of_site(int[:],float[:,:], int)
def get_score_of_site(site, pwm, length):
    score = 0.
    for index in range(length):
        nuc = site[index]
        if nuc < 4:
            score += pwm[nuc][index]
        else:
            score -= 100.
    return score


#pythran export get_scores_of_seq(int[],float[:,:], int, int)
def get_scores_of_seq(seq, pwm, length, number_of_scores):
    scores = np.zeros(number_of_scores * 2, dtype=np.float64)
    for i in range(number_of_scores):
        #forward
        site = seq[i:length + i]
        s = get_score_of_site(site, pwm, length)
        scores[i] = s
        #backward
        if i * -1 == 0:
            site = seq[(length + i) * -1:]
            s = get_score_of_site(site, pwm, length)
            scores[-1] = s
        else:
            site = seq[(length + i) * -1:i * -1]
            s = get_score_of_site(site, pwm, length)
            scores[-i-1] = s
    return scores


#pythran export scaner(int[:,:], float[:,:])
def scaner(promoters, pwm):
    length = pwm.shape[1]
    number_of_promoters = promoters.shape[0]
    number_of_scores = (promoters.shape[1] // 2) - length + 1
    scores = np.zeros((number_of_promoters, number_of_scores * 2), dtype=np.float64)
    for index in range(number_of_promoters):
        seq = promoters[index]
        scores[index] = get_scores_of_seq(seq, pwm, length, number_of_scores)
    return scores


#pythran export calculate_enrichment(float[:,:], float)
#pythran export calculate_enrichment(float[], float)
def calculate_enrichment(scores, threshold):
    counts = np.sum(np.greater_equal(scores, threshold))
    total = scores.shape[0] * scores.shape[1]
    enrichment = counts / total
    if enrichment == 0.0:
        enrichment = 10**(-6)
    return enrichment


#pythran export montecarlo_enrichment(float[:,:], float[:], float[:,:], float[:], float, float)
def montecarlo_enrichment(deg_scores, deg_gc, other_scores, other_gc, threshold, gc_threshold):
    number_of_deg, number_of_scores = deg_scores.shape
    number_of_other = other_scores.shape[0]
    gc_index = np.abs([other_gc - i for i in deg_gc])
    gc_index = gc_index < gc_threshold
    gc_index = [np.where(i)[0] for i in gc_index]
    vec_random_enrichment = np.zeros(1000, dtype=np.float64)
    real_enrichmnet = calculate_enrichment(deg_scores, threshold)
    for i in range(1000):
        index = [random.choice(j) for j in gc_index]
        sample = other_scores[index]
        vec_random_enrichment[i] = calculate_enrichment(sample, threshold)
    random_enrichmnet_mean, random_enrichmnet_std = np.mean(vec_random_enrichment), np.std(vec_random_enrichment)
    z_score = np.abs((real_enrichmnet - random_enrichmnet_mean) / random_enrichmnet_std)
    return z_score, np.log2(real_enrichmnet / random_enrichmnet_mean)


#pythran export montecarlo_fraction(float[:], float[:], float[:], float[:], float, float)
def montecarlo_fraction(deg_scores, deg_gc, other_scores, other_gc, threshold, gc_threshold):
    number_of_deg = deg_scores.shape[0]
    number_of_other = other_scores.shape[0]
    gc_index = np.abs([other_gc - i for i in deg_gc])
    gc_index = gc_index < gc_threshold
    gc_index = [np.where(i)[0] for i in gc_index]
    number_of_deg_with_tfbs = np.sum(np.greater_equal(deg_scores, threshold))
    real_fraction = number_of_deg_with_tfbs / number_of_deg
    vec_random_fraction = np.zeros(1000, dtype=np.float64)
    for i in range(1000):
        index = [random.choice(j) for j in gc_index]
        sample = other_scores[index]
        number_of_other_with_tfbs = np.sum(np.greater_equal(sample, threshold))
        vec_random_fraction[i] = number_of_other_with_tfbs / number_of_deg
    random_fraction_mean, random_fraction_std = np.mean(vec_random_fraction), np.std(vec_random_fraction)
    z_score = np.abs((real_fraction - random_fraction_mean) / random_fraction_std)
    return z_score, np.log2(real_fraction / random_fraction_mean)