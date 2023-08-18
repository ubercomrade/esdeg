import os
import sys
import numpy as np
from operator import itemgetter


## Promoters parser
def complement(seq):
    return seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]


def promoters_parser(path):
    container = []
    gname = ''
    with open(path) as file:
        for line in file:
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
