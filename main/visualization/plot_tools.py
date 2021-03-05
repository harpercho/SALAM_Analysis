import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import matplotlib.pyplot as plt
import sys, os, glob, pickle, struct
import plot_tools

def load_halos_pickle(pickle_path):
    """ Returns pickle as a dictionary"""
    
    data = pickle.load( open( pickle_path , "rb" ))    
    output = dict([(str(k),np.zeros(len(data))) for k in data[1]])
    
    for i in range(len(data)):
        gal_dict = data[i]
        if gal_dict is not None:
            for key, value in output.items():
                try:
                    value[i] = gal_dict[str(key)] + 1e-12
                except:
                    pass
    return output


def plot_mf(a, array, n_bin, zred,vol=60.0,):
    """Plot a mass function. Returns x and y positions"""
    
    hist, bin_edges = np.histogram(np.log10(array),bins=n_bin)

    binmps = np.zeros(len(hist))
    binsize = np.zeros(len(hist))

    for i in np.arange(len(hist)):
        binmps[i] = np.mean([bin_edges[i],bin_edges[i+1]])
        binsize[i] = bin_edges[i+1] - bin_edges[i]
    
    return binmps, hist/(vol**3)/binsize


def filter_list(a,min_lim,max_lim):
    """Remove nans and values outside of range specified"""
    
    def is_valid(elm):
        return elm > min_lim and elm < max_lim and np.isfinite(elm)

    #print(len(b))
    for idx in range(len(a) - 1, -1, -1):
        if not is_valid(a[idx]):
                a = np.delete(a,idx)
    return a


def do_filter(a, b):
    """Filter two lists such that if one entry is invalid, both elements from the lists are removed."""
    
    def is_valid(elm):
        return elm > 1 and np.isfinite(elm)

    #print(len(b))
    for idx in range(len(b) - 1, -1, -1):
        #print(idx)
        if not is_valid(b[idx]):
            a = np.delete(a,idx)
            b = np.delete(b,idx)
    return a,b


def nihao(quantity, z):
    """Load NIHAO pickles"""

    NIHAO_PATH = "/scratch/hc2347/pickles/nihao/pickle_NIHAO.p"

    nihao = pickle.load(open(NIHAO_PATH, 'rb'))
    retval = []
    
    for (k,v) in nihao.items():
        for (zgal, gal) in v.items():
            if zgal < 1:
                zerogal = gal
            if zgal > 3:
                fourgal = gal
        try:
            if z == 0:
                retval.append(zerogal[quantity])
            elif z == 4:
                halo_dict = fourgal
                retval.append(fourgal[quantity])
            else:
                print("Else: " + str(z))
        except:
            retval.append(0)
                
    return np.array(retval)