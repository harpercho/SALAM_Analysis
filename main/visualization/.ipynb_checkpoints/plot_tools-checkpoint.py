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


def load_iso(pickle_path, isostatus = True):    
    data = pickle.load( open( pickle_path , "rb" ))    
    output = dict([(str(k), []) for k in data[1]])
    
    isolated = pickle.load(open("/scratch/hc2347/pickles/iso/00832.p", "rb"))
    
    for gal_dict in data:
        if isolated[gal_dict['ID']] == isostatus:
            for key, value in output.items():
                try:
                    value.append(gal_dict[str(key)])
                except:
                    value.append(0)
    for k, v in output.items():
        output[k] = np.array(v)
    return output


def loadClassify(pickle_path):
    def isInMoster(halo, z):
        y, yerr = moster(halo['mvir'], z)
        return float(y/yerr) <= float(halo['mstar']) <= float(y*yerr)    
    
    def addToOutput(output, gal_dict):
        for key, value in output.items():
            try:
                value.append(gal_dict[str(key)])
            except:
                value.append(0)
                
    data = pickle.load( open( pickle_path , "rb" ))    
    good_output = dict([(str(k), []) for k in data[1]])
    fake_output = dict([(str(k), []) for k in data[1]])
        
    for gal_dict in data:
        if gal_dict is not None:
            if isInMoster(gal_dict, 0.2):
                addToOutput(good_output, gal_dict)
            else:
                gal_dict['mstar'], yerr = moster(gal_dict['mvir'], 0.2)
                addToOutput(fake_output, gal_dict)
  
    for k, v in good_output.items():
        good_output[k] = np.array(v)
    for k, v in fake_output.items():
        fake_output[k] = np.array(v)
        
    return good_output, fake_output


def bin_stats(x, values):
    bin_means, bin_edges, binnumber = stats.binned_statistic(x, values, statistic='mean', bins=25)
    bin_std, bin_edges, binnumber = stats.binned_statistic(x, values, statistic='std', bins=25)
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width/2
    return bin_centers, bin_means, bin_std


def plot_mf(a, array, n_bin, zred,vol=60.0,):
    """Another mass function plotter. This one is less useful but I keep it around in case I need it."""
    
    hist, bin_edges = np.histogram(np.log10(array),bins=n_bin)

    binmps = np.zeros(len(hist))
    binsize = np.zeros(len(hist))

    for i in np.arange(len(hist)):
        binmps[i] = np.mean([bin_edges[i],bin_edges[i+1]])
        binsize[i] = bin_edges[i+1] - bin_edges[i]
    
    return binmps, hist/(vol**3)/binsize


def plot_smf(array, size):
    """Plot mass function. Size refers to volume size."""
    
    mstar = plot_tools.filter_list(array, 1, 1e20)
    logM = np.log10(mstar)
    nbins = 50
    V = size**3
    Phi, edg = np.histogram(logM, bins=nbins)    #Unnormalized histogram and bin edges
    dM = edg[1] - edg[0]
    Max = np.array(edg[0:-1]) + dM/2
    yerr = np.sqrt(Phi)/V/dM
    
    Phi = Phi/V/dM
    return Max, Phi, yerr


def filter_list(a,min_lim,max_lim):
    """Remove nans and values outside of range specified"""
    
    def is_valid(elm):
        return elm > min_lim and elm < max_lim and np.isfinite(elm)

    #print(len(b))
    for idx in range(len(a) - 1, -1, -1):
        if not is_valid(a[idx]):
                a = np.delete(a,idx)
    return a


def do_filter(a, b, lower_limit=1):
    """Filter two lists such that if one entry is invalid, both elements from the lists are removed."""
    
    def is_valid(elm):
        if elm is None:
            return False
        return elm > lower_limit and np.isfinite(elm)

    #print(len(b))
    for idx in range(len(b) - 1, -1, -1):
        #print(idx)
        if not is_valid(b[idx]) or not is_valid(a[idx]) :
            a = np.delete(a,idx)
            b = np.delete(b,idx)
    return a,b


def nihao(quantity, z):
    """Load NIHAO pickles"""

    NIHAO_PATH = "/scratch/hc2347/pickles/nihao/pickle_NIHAO_2.p"

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