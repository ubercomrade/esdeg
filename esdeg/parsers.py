import os
import sys
import lzma
import numpy as np
from operator import itemgetter


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
    

def jaspar_to_pwm(motif):
    pcm = motif.counts
    pcm = dict_to_array(pcm)
    pfm = pcm_to_pfm(pcm)
    pwm = pfm_to_pwm(pfm)
    return pwm

## Promoters parser
def complement(seq):
    return seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]


def promoters_parser(path):
    container = []
    gname = ''
    with lzma.open(path) as file:
        for line in file:
            line = line.decode()
            if not line.startswith('>'):
                seq = line.strip().upper()
                seq += complement(seq)
                container.append((gname, seq))
            else:
                gname = line.strip().split(':')[0][1:]
    container.sort(key=itemgetter(0)) 
    promoters_ids = [i[0] for i in container]
    promoters = [i[1] for i in container]
    return promoters, np.array(promoters_ids)


def fasta_parser(path):
    fasta = []
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                seq = line.strip().upper()
                seq += complement(seq)
                fasta.append(seq)
    return fasta

    
def read_set_of_genes(path):
    container = []
    with open(path) as file:
        for line in file:
            container.append(line.strip())
    return np.array(list(set(container)))
