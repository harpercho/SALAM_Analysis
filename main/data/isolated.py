import pynbody
import numpy as np
import sys
import pickle
from time import perf_counter

# Formally, we only consider haloes that have no other haloes with mass greater than one-fifth of the virial mass within 3 virial radii. 95 per cent of all parent haloes satisfy this constraint, so it does not bias our sample in any significant way. 

def create_coords(catalog):
    '''Creates a reference for position and mass of all halos. Returns as numpy array.'''
    coords = np.zeros((len(catalog), 4))
    print(coords.shape)
    for i, halo in enumerate(catalog):
        coords[i][0], coords[i][1], coords[i][2] = halo.properties['Xc'], halo.properties['Yc'], halo.properties['Zc']
        coords[i][3] = halo.properties['mass']
    return coords
        

def isIsolated(idx):
    ''' Check if halo (by index) is isolated. Returns boolean. '''
    def is_near(other_halo, dist_lim):
        x_, y_, z_ = other_halo[0], other_halo[1], other_halo[2]
        dist = ((x_ - x)**2 + (y_ - y)**2 + (z_ - z)**2)**0.5
        return dist < dist_lim
        
    halo = catalog[idx]
    x, y, z = halo.properties['Xc'], halo.properties['Yc'], halo.properties['Zc']
    dist_lim = 3*halo.properties['Rvir']

    i = 1
    while i < idx:
        other_halo = coords[i - 1]
        if is_near(other_halo, dist_lim):
            print("Found a non-isolated halo. Halo ID", idx)
            return False
        i += 1
    
    # Check for halos farther down
    i =  idx
    while i < idx + 100 and i < len(coords):
        other_halo = coords[i]
        if other_halo[3] > halo.properties['mass']/2:
            if is_near(other_halo, dist_lim):
                print("Found a non-isolated halo. Halo ID", idx)
                return False
        i += 1
    return True

if __name__ == '__main__':

    #path = sys.argv[1]
    path = "/scratch/kld8/simulations/LRZ_Planck60/planck.new.hydro.60_600.00832"
    coordspath = "/scratch/hc2347/references/" + path.split(".")[-1] + "_newcoords.p"
    print("Coordspath is ", coordspath)

    s = pynbody.load(path)
    catalog = s.halos()
    
    try:
        coords = pickle.load(open(coordspath, 'rb'))
    except:
        "Coords pickle does not exist yet."
        coords = create_coords(catalog)
        pickle.dump(coords, open(coordspath, 'wb'))
    
    # idxs = np.random.randint(1, len(catalog), 20)
    # start = perf_counter()
    isolated = {}

    for idx in range(len(catalog)):
        #print("Is halo {} isolated?".format(idx + 1))
        isolated[idx + 1] = isIsolated(idx + 1)
    
    pickle.dump(isolated, open("/scratch/hc2347/pickles/iso/00832_new.p", 'wb'))
    # end = perf_counter())

    # print("The time taken for the loop is {} s".format(end-start))
