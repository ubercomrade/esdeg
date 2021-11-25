import numpy as np
from numba import njit
from operator import itemgetter


### MATRIX PART ###
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
                    line = [float(i) * nsites for i in line.strip().split('\t')]
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
        name, matrix = matrices[index]
        matrix = pcm_to_pfm(matrix)
        matrix = pfm_to_pwm(matrix)
        matrix_length = matrix.shape[1]
        middle_score = to_score(matrix, 0.6)
        container.append([name, matrix, matrix_length, middle_score])
    return container


### FASTA PART ###
def fasta_parser(path):
    container = []
    gname = ''
    seq = ''
    with open(path) as file:
        for line in file:
            if line.startswith('>'):
                if not gname == '':
                    seq += complement(seq)
                    seq = np.array(seq, dtype='c')
                    seq = seq.view(np.uint8)
                    container.append((gname, actg_to_numbers(seq)))
                gname = line.strip().split(':')[0][1:]
                seq = ''
            else:
                seq += line.strip().upper()
    return container


@njit(cache=True)
def actg_to_numbers(seq):
    length = len(seq)
    vec = np.zeros(length, dtype=np.int64)
    for index in range(length):
        if seq[index] == 65:
            vec[index] = 0
        elif seq[index] == 67:
            vec[index] = 1
        elif seq[index] == 71:
            vec[index] = 2
        elif seq[index] == 84:
            vec[index] = 3
        else:
            vec[index] = 4
    return vec
    
    
def complement(seq):
    return seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]