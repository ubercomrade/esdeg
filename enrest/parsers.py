import sys
import lzma
import numpy as np
from operator import itemgetter
from enrest.speedup import work_with_seq

def read_list_of_matrix_hocomoco(path):
    container = []
    matrix = []
    with open(path) as file:
        matrix_name = file.readline().strip()[1:]
        for line in file:
            if line.startswith(">"):
                container.append((matrix_name, np.array(matrix).T))
                matrix_name = line.strip()[1:]
                matrix = []
            else:
                line = [float(i) for i in line.strip().split('\t')]
                matrix.append(line)
    container.append((matrix_name, np.array(matrix).T))
    return container


def read_list_of_matrix_meme(path):
    container = []
    matrix = []
    with open(path) as file:
        for line in file:
            if line.startswith("MOTIF"):
                matrix_name = line.split()[1]
                for line in file:
                    if line.strip() != "":
                        length, nsites = map(int, itemgetter(*[5,7])(line.strip().split()))
                        break
                for index, line in enumerate(file):
                    line = [float(i) * nsites for i in line.strip().split()]
                    matrix.append(line)
                    if index == length - 1:
                        container.append((matrix_name, np.array(matrix).T))
                        matrix = []
                        break
    return container


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


def matrices_parser(path, f="meme"):
    container = []
    try:
        if f == "meme":
            matrices = read_list_of_matrix_meme(path)
        elif f == "hocomoco":
            matrices = read_list_of_matrix_hocomoco(path)
        else:
            sys.exit("Wrong format file")
    except:
        sys.exit("Can't read file with matrices")
    for index in range(len(matrices)):
        name, pcm = matrices[index]
        pfm = pcm_to_pfm(pcm)
        pwm = pfm_to_pwm(pfm)
        matrix_length = pwm.shape[1]
        container.append([name, pwm, pfm, matrix_length])
    return container


def promoters_parser(path):
    promoters = []
    promoters_ids = []
    gname = ''
    with lzma.open(path) as file:
        for line in file:
            line = line.decode()
            if not line.startswith('>'):
                seq = work_with_seq(line)
                promoters.append(seq)
                promoters_ids.append(gname)
            else:
                gname = line.strip().split(':')[0][1:]
    return np.array(promoters), promoters_ids


def fasta_parser(path):
    container = []
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                seq = work_with_seq(line)
                container.append(seq)
    return np.array(container)

    
def read_set_of_genes(path):
    container = []
    with open(path) as file:
        for line in file:
            container.append(line.strip())
    return set(container)
