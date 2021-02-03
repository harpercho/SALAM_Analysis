# A module for functions to return required quantities.
# Halo can be a halo copy. Storage must be dict.
    

def getMasses(halo, storage):
    mvir = np.sum(halo['mass'].in_units('Msol'))
    mstar = np.sum(halo.star['mass'].in_units('Msol'))
    mgas = np.sum(halo.gas['mass'].in_units('Msol'))
    mdm = np.sum(halo.dm['mass'].in_units('Msol'))
    
    storage.update({'mstar':mstar, 'mgas':mgas})
    
    #print("Obtained masses", flush=True)
    
def getColdGas(halo, storage, temp):
    # Halo can be a halo copy. Storage must be dict.
    coolgasf = filt.And(filt.LowPass('temp',temp),filt.HighPass('rho','0.03 m_p cm^-3'))
    
    storage['mgascool'] = np.sum(halo.gas[coolgasf]['mass'].in_units('Msol'))
    
    #print("Obtained coldgas", flush = True)
    
def getParticleInfo(halo, storage):
    
    npart = len(halo)
    nstar = len(halo.star)
    ngas = len(halo.gas)
    ndm = len(halo.dm)
    
    storage.update({'npart': npart, 'nstar': nstar, 'ngas': ngas, 'ndm':ndm})
    
    #print("Obtained particle info", flush = True)

def getSFR(halo, storage, Myr):
    # Halo can be a copy
    # Myr must be an int or float
    fifmyrf = filt.LowPass('age', str(Myr) + ' Myr')
    storage['sfr_'+str(Myr)] = np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / (Myr*10**6)
    
    #print("Obtained SFR", flush = True)
    
def getMetallicity(halo, storage):
    mgas = np.sum(halo.gas[coolgasf]['mass'].in_units('Msol'))
    
    if mgas > 0:
        zgas = np.sum(halo.g[coolgasf]['mass'].in_units('Msol')*halo.g[coolgasf]['metals'])/mgas
    else:
        zgas = 0
    storage['z_gas'] = zgas
    
    #print("Obtained metallicity", flush = True)
    
def getStellarMetallicity(halo, storage, mstar):
    if mstar > 0:
        zstar = np.sum(halo.s['mass'].in_units('Msol')*halo.s['metals'])/mstar
    else:
        zstar = 0
    storage['z_star'] = zstar
    #print("Obtained STELLAR METALLICITY", flush = True)

def getOxygenAbundance(halo, storage, temp):
    '''
    According to Tremonti et al 2004
    '''
    
    coolgasf = filt.And(filt.LowPass('temp',temp),filt.HighPass('rho','0.03 m_p cm^-3'))
    
    try:
        storage['oxh'] = np.log10(np.sum(halo.g[coolgasf]['OxMassFrac'])/(16*np.sum(halo.g[coolgasf]['hydrogen']))) + 12
    except:
        print("Unable to obtain oxh")
        return None