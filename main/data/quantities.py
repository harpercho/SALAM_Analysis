# A module for functions to return required quantities.
# Halo can be a halo copy. Storage must be dict.
    
import os
import sys
import numpy as np
import pickle
import pynbody
import pynbody.filt as filt
import multiprocessing
import concurrent.futures

def getMasses(halo, storage):
    mvir = np.sum(halo['mass'].in_units('Msol'))
    mstar = np.sum(halo.star['mass'].in_units('Msol'))
    mgas = np.sum(halo.gas['mass'].in_units('Msol'))
    mdm = np.sum(halo.dm['mass'].in_units('Msol'))
    
    storage.update({'mstar':mstar, 'mgas':mgas})
    
    
def getColdGas(halo, storage):
    # Halo can be a halo copy. Storage must be dict.
    temp = 1.5e4
    coolgasf = filt.And(filt.LowPass('temp',temp),filt.HighPass('rho','0.03 m_p cm^-3'))
    
    storage['mgascool'] = np.sum(halo.gas[coolgasf]['mass'].in_units('Msol'))
    
    
def getParticleInfo(halo, storage):
    
    npart = len(halo)
    nstar = len(halo.star)
    ngas = len(halo.gas)
    ndm = len(halo.dm)
    
    if nstar == 0:
        print("This halo has no stars.")
    
    storage.update({'npart': npart, 'nstar': nstar, 'ngas': ngas, 'ndm':ndm})
    

def getSFR(halo, storage, Myr):
    
    fifmyrf = filt.LowPass('age', str(Myr) + ' Myr')
    storage['sfr_'+str(Myr)] = np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / (Myr*10**6)
    
    
def getMetallicity(halo, storage):
    mgas = np.sum(halo.gas[coolgasf]['mass'].in_units('Msol'))
    
    if mgas > 0:
        zgas = np.sum(halo.g[coolgasf]['mass'].in_units('Msol')*halo.g[coolgasf]['metals'])/mgas
    else:
        zgas = 0
    storage['z_gas'] = zgas
    
    #print("Obtained metallicity", flush = True)
    
def getStellarMetallicity(halo, storage):
    mstar = np.sum(halo.star['mass'].in_units('Msol'))
    if mstar > 0:
        zstar = np.sum(halo.s['mass'].in_units('Msol')*halo.s['metals'])/mstar
    else:
        zstar = 0
    storage['z_star'] = zstar

def getOxygenAbundance(halo, storage):
    '''
    According to Tremonti et al 2004
    '''
    
    temp = 1.5e4
    coolgasf = filt.And(filt.LowPass('temp',temp),filt.HighPass('rho','0.03 m_p cm^-3'))
    
    try:
        coolgash = np.sum(halo.g[coolgasf]['hydrogen'])
        if coolgash == 0:
            oxh = 0
            print("No cool gas for halo with stellar radius " + str(np.max(halo.star['r'].in_units('kpc'))))
        else:
            oxh = np.log10(np.sum(halo.g[coolgasf]['OxMassFrac'])/(16*coolgash)) + 12
        storage['oxh'] = oxh
    except Exception as e:
        print(e)
        print("Unable to obtain oxh")
        return None