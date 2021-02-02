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

def getMasses(halo, storage):
    # Halo can be a halo copy. Storage must be dict.
    #mvir = np.sum(halo['mass'].in_units('Msol'))
    mstar = np.sum(halo.star['mass'].in_units('Msol'))
    mgas = np.sum(halo.gas['mass'].in_units('Msol'))
    #mdm = np.sum(halo.dm['mass'].in_units('Msol'))
    
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
        
def make_future_runner(idx):
    
    # The script run by worker processes.
    print("{}: Working on halo {}/{}".format(multiprocessing.current_process().name, idx ,num_halos), flush = True)
    
    # If there is no halo dictionary, create a halo dictionary and give it an ID property.
    halo_dict = {}
    
    # Load a copy of the target halo.
    halo = halo_catalogue.load_copy(idx)

    # Center halo
    center(halo) #applies a transformation
    smallhalo = radius_cut(halo, idx)
    
    # Give physical units.
    halo.physical_units()
    
    coldgastemp = 1.5e4
    
    # From quantities.py. This file contains all the actual quantity obtaining utilities.
    # Each of these functions require a storage dictionary to update/add into.
    print("starting obtaining process for halo " + str(idx),flush=True)
    getParticleInfo(smallhalo, halo_dict)
    getMasses(halo, halo_dict)
    getSFR(smallhalo, halo_dict, 10)
    getSFR(smallhalo, halo_dict, 100)
    getColdGas(halo, halo_dict, coldgastemp)
    getOxygenAbundance(halo, halo_dict, coldgastemp)
    getStellarMetallicity(halo, halo_dict, halo_dict['mstar'])
    
    
    halo_dict['ID'] = idx
    halo_dict['mvir'] = np.sum(halo['mass'].in_units('Msol'))
    halo_dict['mdm'] = np.sum(halo.dm['mass'].in_units('Msol'))
    
    print('Halo ' + str(idx) + ' complete')
    
    return halo_dict

def main(path, step):
    
    PROC_N = 4 # number of worker processes
    
    # Check if there are halos in the snapshot
    if num_halos <= 1.5:
        print("No halos in this snapshot.")
        return None
    
    # Check if a filtered index list already exists. If not, make a new one.
    if os.path.exists(LIST_PATH):
        print('Found existing filtered halo list.')
        halolist = pickle.load(open(LIST_PATH,'rb'))
    else:
        print('No halo list found. Creating new list.')
        halolist = makeFilteredHalos(halo_catalogue)
        pickle.dump(halolist, open(LIST_PATH,'wb'))
            
    args = [] # list of tuple arguments to be passed to the executor map function
    for halo_idx in halolist:
            args.append([halo_idx, {}])
    
    # Check if there is an existing output file and create a full list of halodict storages
    if os.path.exists(OUTPUT):
        outd = pickle.load(open(OUTPUT,'rb'))
        for item in outd:
            try:
                i = halolist.index(item['ID'])
                args[i][1] = item
            except ValueError:
                pass

    # Distribute
    with concurrent.futures.ProcessPoolExecutor(max_workers=PROC_N) as executor:
        #result = executor.map(make_future_runner, range(1,15))
        result = executor.map(make_future_runner, halolist)
        flat_results = list(result)
    
    # Save it all
    pickle_path = "/scratch/hc2347/pickles/{}".format(box)
    if os.path.isdir(pickle_path) is False:
        os.mkdir(pickle_path)
    
    pickle.dump( flat_results, open( OUTPUT, 'wb' ))
    print('Done')

if __name__ == '__main__':
    
    '''
    Usage: python main/data/preload.py /scratch/kld8/simulations/LRZ_Planck60 planck.new.hydro.60_600.00832
    '''
    
    path = sys.argv[1]
    step = sys.argv[2]
    
    box = step.split('.')[3].split('_')[0]
    
    # The path for the sim we want to analyze.
    SIM_FILE = os.path.join(path, step)
    LIST_PATH = '/scratch/hc2347/references/' + box + '_' + step.split('.')[-1] + '_filtered_halos.p'


    pynbody.config['threading'] = False
    pynbody.config['number_of_threads'] = 1
    print("Main threads: " + str(pynbody.config['number_of_threads']), flush=True)
    
    sim = pynbody.load(SIM_FILE)
    halo_catalogue = sim.halos()
    halo_properties = sim.halos(dummy=True)

    zred = sim.properties['z']    
    num_halos = len(halo_catalogue)
    
    OUTPUT = os.path.join("/scratch/hc2347/pickles/{}".format(box),'centered_z{:.3f}.p'.format(zred))
    
    
    main(path, step)
