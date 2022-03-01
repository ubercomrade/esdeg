import numpy as np
from numba import njit, float64, int64


@njit(float64(int64[:],float64[:,:], int64),fastmath=True, cache=True)
def get_score_of_site(site, pwm, length):
    score = 0.
    for index in range(length):
        nuc = site[index]
        if nuc < 4:
            score += pwm[nuc][index]
        else:
            score -= 100.
    return score


@njit(float64[:](int64[:],float64[:,:], int64, int64), cache=True)
def get_scores_of_seq(seq, pwm, length, number_of_scores):
    scores = np.zeros(number_of_scores * 2)
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


def scaner(promoters, pwm):
    container = []
    length = pwm.shape[1]
    number_of_scores = int(promoters[0][1].shape[0] / 2 - length + 1)
    for gname, seq in promoters:
        scores = get_scores_of_seq(seq, pwm, length, number_of_scores)
        container.append((gname, scores))
    return container