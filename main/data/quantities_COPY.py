# A module for functions to return required quantities.
    
import os
import sys
import numpy as np
import pickle
import pynbody
import pynbody.filt as filt
import multiprocessing
import concurrent.futures
from haloprep import *

def getMasses(halo):
    mvir = np.sum(halo['mass'].in_units('Msol'))
    mstar = np.sum(halo.star['mass'].in_units('Msol'))
    mgas = np.sum(halo.gas['mass'].in_units('Msol'))
    mdm = np.sum(halo.dm['mass'].in_units('Msol'))
    
    return mvir, mstar, mgas, mdm
    
    
def getColdGas(halo, temp):
    """
    Based on Peeples 2014
    Set smallify to False if the halo being passed is already cut according to radius.
    Temperature sets the lowpass limit.
    """
    coolgasf = filt.And(filt.LowPass('temp',temp),filt.HighPass('rho','0.03 m_p cm^-3'))
    mgascool = np.sum(halo.gas[coolgasf]['mass'].in_units('Msol'))
    
    return mgascool
    
def getParticleInfo(halo):
    npart = len(halo)
    nstar = len(halo.star)
    ngas = len(halo.gas)
    ndm = len(halo.dm)
    
    if nstar == 0:
        print("This halo has no stars.")

    print("npart: {}, nstar: {}, ngas: {}, ndm: {}".format(npart, nstar, ngas, ndm))
    
    return npart, nstar, ngas, ndm
    

def getSFR(halo, Myr):
    fifmyrf = filt.LowPass('age', str(Myr) + ' Myr')
    return np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / (Myr*10**6)
    
    
def getStellarMetallicity(halo, mstar):
    if mstar > 0:
        zstar = np.sum(halo.s['mass'].in_units('Msol')*halo.s['metals'])/mstar
    else:
        zstar = 0
    return zstar

    
def getOxygenAbundance(halo, smallify = False):
    '''
    According to Tremonti et al 2004.
    Set smallify to false if the halo is already cut.
    '''
    
    temp = 1.5e4
    coolgasf = filt.And(filt.LowPass('temp',temp),filt.HighPass('rho','0.03 m_p cm^-3'))
    
    # Tremonti 
    if smallify:
        try:
            halo = half_stellar_radius(halo)
        except Exception as e:
            print(e)

    try:
        coolgash = np.sum(halo.g[coolgasf]['hydrogen'])
        if coolgash == 0:
            oxh = 0
            print("No cool gas for halo with stellar radius " + str(np.max(halo.star['r'].in_units('kpc'))))
        else:
            oxh = np.log10(np.sum(halo.g[coolgasf]['OxMassFrac'])/(16*coolgash)) + 12
        return oxh
    except Exception as e:
        print(e)
        print("Unable to obtain oxh")
        return None
