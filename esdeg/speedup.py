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


#pythran export count_sites(float[:,:], float[:,:])
def count_sites(scores, threshold_table):
    counts = np.zeros((scores.shape[0], threshold_table.shape[0]), dtype=np.int64)
    for i in range(threshold_table.shape[0]):
        counts[:,i] = np.sum(np.greater_equal(scores, threshold_table[i][0]), axis=1)
    return counts