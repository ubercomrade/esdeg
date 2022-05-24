# $ pythran -DUSE_XSIMD -fopenmp -march=native speedup.py

import numpy as np
#pythran export sort(float[])
def sort(scores):
    return np.sort(scores)[::-1]
    

#pythran export complement(str)
def complement(seq):
    seq = list(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper())
    reversed_seq = ''.join(reversed(seq))
    return reversed_seq


#pythran export seq_to_int(str)
def seq_to_int(seq):
    length = len(seq)
    vec = np.zeros(length, dtype=np.int64)
    for index in range(length):
        if seq[index] == 'A':
            vec[index] = 0
        elif seq[index] == 'C':
            vec[index] = 1
        elif seq[index] == 'G':
            vec[index] = 2
        elif seq[index] == 'T':
            vec[index] = 3
        else:
            vec[index] = 4
    return vec
 
    
#pythran export work_with_seq(str)
def work_with_seq(seq):
    seq = seq.strip().upper()
    seq = complement(seq)
    seq = seq_to_int(seq)
    return seq
    

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
    # omp parallel for
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


#pythran export montecarlo_enrichment(float[:,:], float[:,:], float)
def montecarlo_enrichment(deg_scores, other_scores, threshold):
    number_of_deg = len(deg_scores)
    number_of_other = len(other_scores)
    #indexes = np.arange(number_of_other)
    vec_random_enrichment = np.zeros(1000, dtype=np.float64)
    real_enrichmnet = calculate_enrichment(deg_scores, threshold)
    for i in range(1000):
        #sample_indexes = np.random.choice(indexes, number_of_deg)
        sample_indexes = np.random.randint(0, number_of_other, number_of_deg)
        sample = other_scores[sample_indexes]
        vec_random_enrichment[i] = calculate_enrichment(sample, threshold)
    random_enrichmnet_mean, random_enrichmnet_std = np.mean(vec_random_enrichment), np.std(vec_random_enrichment)
    z_score = abs((real_enrichmnet - random_enrichmnet_mean) / random_enrichmnet_std)
    return z_score, np.log2(real_enrichmnet / random_enrichmnet_mean)


#pythran export montecarlo_fraction(float[], float[], float)
def montecarlo_fraction(deg_scores, other_scores, threshold):
    number_of_deg = len(deg_scores)
    number_of_other = len(other_scores)
    number_of_deg_with_tfbs = np.sum(np.greater_equal(deg_scores, threshold))
    real_fraction = number_of_deg_with_tfbs / number_of_deg
    #indexes = np.arange(number_of_other)
    vec_random_fraction = np.zeros(1000, dtype=np.float64)
    # omp parallel for
    for i in range(1000):
        #sample_indexes = np.random.choice(indexes, number_of_deg)
        sample_indexes = np.random.randint(0, number_of_other, number_of_deg)
        sample = other_scores[sample_indexes]
        number_of_other_with_tfbs = np.sum(np.greater_equal(sample, threshold))
        vec_random_fraction[i] = number_of_other_with_tfbs / number_of_deg
    random_fraction_mean, random_fraction_std = np.mean(vec_random_fraction), np.std(vec_random_fraction)
    z_score = abs((real_fraction - random_fraction_mean) / random_fraction_std)
    return z_score, np.log2(real_fraction / random_fraction_mean)
