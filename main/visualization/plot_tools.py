import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
import matplotlib.pyplot as plt
import sys, os, glob, pickle

def load_halos_pickle(pickle_path):
    
    '''
    Returns pickle as dictionary
    '''
        
    data = pickle.load( open( pickle_path , "rb" ))
    
    output = dict([(str(k),np.zeros(len(data))) for k in data[0]])
    
    for i in range(len(data)):
        gal_dict = data[i]
        
        if gal_dict is not None:
            for key, value in output.items():
                try:
                    value[i] = gal_dict[str(key)]
                except:
                    pass
    return output

def plot_mf(a, array, n_bin, zred,vol=60.0,):
    hist, bin_edges = np.histogram(np.log10(array),bins=n_bin)

    binmps = np.zeros(len(hist))
    binsize = np.zeros(len(hist))

    for i in np.arange(len(hist)):
        binmps[i] = np.mean([bin_edges[i],bin_edges[i+1]])
        binsize[i] = bin_edges[i+1] - bin_edges[i]
    a.set_ylim(10**-6,1)
    return a.plot(binmps,hist/(vol**3)/binsize,'.',label=str(zred),alpha=0.8)

def do_filter(a, b):
    
    def is_valid(elm):
        return elm > 1 and np.isfinite(elm)

    #print(len(b))
    for idx in range(len(b) - 1, -1, -1):
        #print(idx)
        if not is_valid(b[idx]):
                a = np.delete(a,idx)
                b = np.delete(b,idx)
    return a,b

def filter_list(data_array):
    res = []
    counter = 0
    for i in range(len(data_array)):
        if data_array[i] > 1.5:
            res.append(data_array[i])
        else:
            counter += 1
    print("Removed " + str(counter) + " from array.")
    return np.asarray(res)

