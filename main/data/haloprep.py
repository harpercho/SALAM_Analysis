import os
import sys
import numpy as np
import pickle
import pynbody
import pynbody.filt as filt
import multiprocessing
import concurrent.futures

def isHostHalo(h):
    return h.properties['hostHalo'] <= 0
    
def isBigHalo(h,cut):
    return len(h) > cut
    
def makeFilteredHalos(h):
    
    '''
    Creates a index list of halos that pass our filter.
    Default is 100 particles (generic) and must be hosthalo.
    
    If WriteFile is true, it writes this into a pickle file in ../references.
    '''
    
    # Set minimum number of particles
    P_MIN = 100
    
    halolist = []
    failed = 0
    
    for i,halo in enumerate(h):
        if isHostHalo(halo) and isBigHalo(halo,P_MIN):
            halolist.append(i+1)
            print("Appended " + str(i))
        else:
            print("Removed " + str(i))
            failed +=1
            
    print('Removed ' + str(failed) + 'halos')
    
    return halolist

def center(halo):
    try:
        # Make vel=True only if velocity is useful in the quantities being obtained.
        # By default uses 'ssc' mode according to pynbody.config
        centering = pynbody.analysis.halo.center(halo, vel=False)
    except Exception as e:
        print("Could not center: {}".format(e), flush=True)
    
def radius_cut(halo, idx):
    try:
        RVIR = pynbody.array.SimArray(np.max(halo['r'].in_units('kpc')),'kpc')
    except Exception as e:
        print("No RVIR. Trying AHF method.")
        print(e)
        try:
            RVIR = halo_properties[idx].properties['Rvir']
        except:
                print("No RVIR.")
                return halo
                
    diskf = filt.Sphere(str(RVIR*0.2) +' kpc') #20% Virial Radius
    
    return halo[diskf]

def half_stellar_radius(halo):
    half_sm = np.sum(halo.star['mass'].in_units('Msol'))* 0.5

    max_high_r = np.max(halo.star['r'].in_units('kpc'))
    #print(max_high_r)
    test_r = 0.5 * max_high_r
    testrf = filt.LowPass('r', test_r)
    min_low_r = 0.0
    test_sm = np.sum(halo[testrf].star['mass'].in_units('Msol'))
    it = 0
    while ((np.abs(test_sm - half_sm) / half_sm) > 0.01):
        it = it + 1
        if (it > 20):
            break

        if (test_sm > half_sm):
            test_r = 0.5 * (min_low_r + test_r)
        else:
            test_r = (test_r + max_high_r) * 0.5
        testrf = filt.LowPass('r', test_r)
        test_sm = np.sum(halo[testrf].star['mass'].in_units('Msol'))

        if (test_sm > half_sm):
            max_high_r = test_r
        else:
            min_low_r = test_r

    #print("Half stellar radius found as: {}".format(test_r))
    diskf = filt.Sphere(str(test_r) +' kpc')
    
    return halo[diskf]