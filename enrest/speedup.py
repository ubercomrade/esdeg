import numpy as np


#pythran export sort(float[])
def sort(scores):
    return np.sort(scores)[::-1]
