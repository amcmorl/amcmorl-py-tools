import numpy as n
from scipy.ndimage import label

def biggest(bivol):
    '''return only largest connected structure in the binary volume
    bivol'''
    lbls = label(bivol)
    sizes = n.asarray(n.histogram(lbls[0], bins=range(1,lbls[1]+1)))
    maxidx = sizes[1][sizes[0] == sizes[0].max()].squeeze()
    return lbls[0] == maxidx
