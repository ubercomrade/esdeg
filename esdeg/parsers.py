import os
import sys
import numpy as np
from operator import itemgetter


def promoters_parser(path):
    container = []
    gname = ''
    letters = {'A', 'C', 'G', 'T'}
    with open(path) as file:
        for line in file:
            if line.startswith('>'):
                if gname != '':
                    container.append((gname, seq))
                gname = line.strip().split(':')[0][1:]
                seq = ''
            else:
                seq += ''.join([l if l in letters else 'N' for l in line.strip().upper()])
        container.append((gname, seq))
    container.sort(key=itemgetter(0))
    promoters_ids = [i[0] for i in container]
    promoters = [i[1] for i in container]
    return promoters, np.array(promoters_ids)


def read_set_of_genes(path):
    container = []
    with open(path) as file:
        for line in file:
            container.append(line.strip())
    return np.array(list(set(container)))

# def complement(seq):
#     return seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]

# def fasta_parser(path):
#     fasta = []
#     with open(path) as file:
#         for line in file:
#             if not line.startswith('>'):
#                 seq = line.strip().upper()
#                 seq += complement(seq)
#                 fasta.append(seq)
#     return fasta

#OLD
# def promoters_parser(path):
#     container = []
#     gname = ''
#     with open(path) as file:
#         for line in file:
#             if not line.startswith('>'):
#                 seq = line.strip().upper()
#                 seq += complement(seq)
#                 container.append((gname, seq))
#             else:
#                 gname = line.strip().split(':')[0][1:]
#     container.sort(key=itemgetter(0))
#     promoters_ids = [i[0] for i in container]
#     promoters = [i[1] for i in container]
#     return promoters, np.array(promoters_ids)


#OLD
# def promoters_parser(path):
#     container = []
#     gname = ''
#     with open(path) as file:
#         for line in file:
#             if line.startswith('>'):
#                 if gname != '':
#                     container.append((gname, seq + complement(seq)))
#                 gname = line.strip().split(':')[0][1:]
#                 seq = ''
#             else:
#                 seq += line.strip().upper()
#         container.append((gname, seq + complement(seq)))
#     container.sort(key=itemgetter(0))
#     promoters_ids = [i[0] for i in container]
#     promoters = [i[1] for i in container]
#     return promoters, np.array(promoters_ids)
