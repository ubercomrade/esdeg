import json
import numpy as np
from importlib import resources
from pyjaspar import jaspardb
from operator import itemgetter

#read promoters
def promoters_parser(path):
    container = []
    gname = ''
    seq = ''
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


#Read motif DB
def dict_to_array(motif):
    motif = [motif[i] for i in motif.keys()]
    return np.array(motif)


def pfm_to_pwm(pfm):
    background = 0.25
    pwm = np.log(pfm / background)
    return pwm


def pcm_to_pfm(pcm):
    number_of_sites = pcm.sum(axis=0)
    nuc_pseudo = 0.25
    pfm = (pcm + nuc_pseudo) / (number_of_sites + 1)
    return pfm


def check_motif_annotaion(ann):
    if len(ann) == 0:
        ann = 'NA'
    elif len(ann) == 1:
        ann = ann[0]
    else:
        ann = '::'.join(ann)
    return ann


def read_motifs_from_db(motif_db, taxon):
    print(f'Read motifs DB: {motif_db}')
    container = []
    if motif_db == 'hocomoco':
        hocomoco_path = resources.files('esdeg').joinpath('hocomoco/H12CORE_annotation.jsonl')
        with open(hocomoco_path) as file:
            motifs = file.readlines()
        motifs = [json.loads(i) for i in motifs]
        motifs = [i for i in motifs if i['length'] >= 6]
        number_of_motifs = len(motifs)
        for motif_data in motifs:
            motif_id = motif_data['name']
            if 'HUMAN' in motif_data['masterlist_info']['species']:
                tf_name = motif_data['masterlist_info']['species']['HUMAN']['gene_symbol'].upper()
            else:
                tf_name = motif_data['masterlist_info']['species']['MOUSE']['gene_symbol'].upper()
            tf_class = motif_data['masterlist_info']['tfclass_class']
            tf_family = motif_data['masterlist_info']['tfclass_family']

            pcm = np.array(motif_data['pcm']).T
            pfm = pcm_to_pfm(pcm)
            pfm = np.concatenate((pfm, np.min(pfm, axis=0).reshape(1, pfm.shape[1])), axis=0)
            pwm = pfm_to_pwm(pfm)
            pwm = pwm.astype(np.float64)

            container.append(
                {'motif_id': motif_id,
                'tf_name': tf_name,
                'tf_class': tf_class,
                'tf_family': tf_family,
                'pfm': pfm,
                'pwm': pwm,
                'length': pfm.shape[1]}
            )


    elif motif_db == 'jaspar':
        #plants vertebrates insects urochordates nematodes fungi
        jdb_obj = jaspardb(release='JASPAR2024')
        motifs = jdb_obj.fetch_motifs(collection = 'CORE',tax_group = [taxon], min_length=6)
        number_of_motifs = len(motifs)
        for motif_data in motifs:
            motif_id = motif_data.matrix_id
            tf_name = motif_data.name
            tf_class = check_motif_annotaion(motif_data.tf_class)
            tf_family = check_motif_annotaion(motif_data.tf_family)

            pcm = motif_data.counts
            pcm = dict_to_array(pcm)
            pfm = pcm_to_pfm(pcm)
            pfm = np.concatenate((pfm, np.min(pfm, axis=0).reshape(1, pfm.shape[1])), axis=0)
            pwm = pfm_to_pwm(pfm)
            pwm = pwm.astype(np.float64)

            container.append(
                {'motif_id': motif_id,
                'tf_name': tf_name,
                'tf_class': tf_class,
                'tf_family': tf_family,
                'pfm': pfm,
                'pwm': pwm,
                'length': pfm.shape[1]}
            )

    print(f'Number of matrices = {number_of_motifs}')
    print('-'*30)
    return container
