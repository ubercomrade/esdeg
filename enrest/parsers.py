import numpy as np
import logomaker
import matplotlib.pyplot as plt
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
        name, pcm = matrices[index]
        pfm = pcm_to_pfm(pcm)
        pwm = pfm_to_pwm(pfm)
        matrix_length = pwm.shape[1]
        middle_score = to_score(pwm, 0.6)
        container.append([name, pwm, pfm, matrix_length, middle_score])
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


### SET GENES ###
def read_set_of_genes(path):
    container = []
    with open(path) as file:
        for line in file:
            container.append(line.strip())
    return set(container)


### PLOTS ###
def plot_bs_distribution(set_scores, threshold_table, write_path, length=2000, window=20):
    fig, ax = plt.subplots(figsize=(4, 2), dpi=200)
    number_of_genes = len(set_scores)
    for index, level in zip(range(0, len(threshold_table)), ['LOW', 'MIDDLE', 'HIGH']):
        threshold, fpr = threshold_table[index]
        distribution = np.sum(np.greater_equal(set_scores, threshold), axis=0)
        add_zeros = int(length - len(distribution) // 2)
        distribution = np.concatenate([distribution.reshape(2, len(distribution) // 2),
                                       np.zeros(add_zeros * 2, dtype=np.int64).reshape(2, add_zeros)],
                                      axis=1)
        distribution[1] = distribution[1][::-1]
        distribution = np.sum(distribution, axis=0)
        distribution = np.sum(distribution.reshape(window, length // window), axis=0)
        ax.plot(np.arange(len(distribution)), distribution / number_of_genes, label=level)
        ax.set_xlim([0,len(distribution) - 1])
        ax.axes.xaxis.set_visible(False)
        ax.legend()
    plt.savefig(write_path, format="jpg", dpi=200, bbox_inches='tight', pad_inches = 0.1)
    plt.close()
    pass


def plot_logo(pfm, write_path):
    bits = np.sum(pfm * np.log2(pfm / 0.25), axis=0)
    bits = bits * pfm
    bits_df = pd.DataFrame(bits.T)
    bits_df.columns = ['A', 'C', 'G', 'T']
    fig, ax = plt.subplots(figsize=(10, 2.5), dpi=200)
    ss_logo = logomaker.Logo(bits_df,
                             width=.9,
                             vpad=.02,
                             stack_order='big_on_top',
                             font_name='DejaVu Sans',
                             show_spines=False,
                             ax=ax)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.savefig(write_path, format="jpg", dpi=200, bbox_inches='tight', pad_inches = 0.1)
    plt.close()
    pass