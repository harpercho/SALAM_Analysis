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
        if isHostHalo(halo[i]) and isBigHalo(halo[i],P_MIN):
            halolist.append(i+1)
        else:
            failed +=1
            
    print('Removed ' + str(failed) + 'halos')
    
    return halolist
  
def make_future_runner(idx):
    print("{}: Working on halo {}/{}".format(multiprocessing.current_process().name, idx ,num_halos), flush = True)
    
    halo_dict = {}
    print(halo_catalogue[idx])
    print("Loading copy for halo " + str(idx), flush = True)
    halo = halo_catalogue.load_copy(idx)
    print("Copy successfully loaded.", flush = True)
    
#     try:
#         RVIR = halo_catalogue[idx].properties['Rvir']
#     except:
#         print('no RVIR',flush=True)
#         return None
    
#     print("finding center of " + str(idx), flush=True)
#     diskf = filt.Sphere(str(RVIR*0.2) +' kpc',pynbody.analysis.halo.center_of_mass(halo)) #20% Virial Radius
#     print("center of " + str(idx) + "found", flush = True)
    
#     halo = halo[diskf]
    
    halo.physical_units()
    
    #getParticleInfo(halo, halo_dict)
    npart = len(halo)
    nstar = len(halo.star)
    ngas = len(halo.gas)
    ndm = len(halo.dm)
    
    mvir = np.sum(halo['mass'].in_units('Msol'))
    mstar = np.sum(halo.star['mass'].in_units('Msol'))
    mgas = np.sum(halo.gas['mass'].in_units('Msol'))
    mdm = np.sum(halo.dm['mass'].in_units('Msol'))
    
    fifmyrf = filt.LowPass('age', str(10) + 'Myr')
    sfr10 = np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / (10*10**6)
    
    fifmyrf = filt.LowPass('age', str(100) + 'Myr')
    sfr100 = np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / (100*10**6)
    
    if mgas > 0:
        zgas = np.sum(halo.g['mass'].in_units('Msol')*halo.g['metals'])/mgas
    else:
        zgas = 0
    
    coolgasf = filt.And(filt.LowPass('temp',temp),filt.HighPass('rho','0.03 m_p cm^-3'))
    mgascool = np.sum(halo.gas[coolgasf]['mass'].in_units('Msol'))
    
    try:
        gas_part = halo.g['mass'].in_units('Msol')
        ox_tot = halo.g['OxMassFrac']*gas_part
        h_tot = halo.g['hydrogen']*gas_part
        
        if h_tot == 0:
            oxygen_abundance = 0
        else:
            oxygen_abundance = np.log10(np.sum(ox_tot/16/h_tot)) + 12
       # print('Oxygen_abundance {}'.format(oxygen_abundance))
    except:
        print('Missing Oxygen Abundance for halo {}'.format(idx))
        oxygen_abundance = 0
    

    if mstar > 0:
        zstar = np.sum(halo.s['mass'].in_units('Msol')*halo.s['metals'])/mstar
    else:
        zstar = 0
    
    halo_dict = {
        'npart':npart,
        'nstar':nstar,
        'ngas':ngas,
        'ndm':ndm,
        'mvir': mvir,
        'mstar':mstar,
        'mgas': mgas,
        'mdm': mdm,
        'sfr10': sfr10,
        'sfr100': sfr100,
        'z_gas': zgas,
        'mgascool': mgascool,
        'oxh': oxygen_abundance,
        'z_star': zstar
        
    }
    
    return halo_dict

def main(path, step, pickle_path=None):
    
    PROC_N = 20  # number of worker processes
        
    print("Loading sim {}".format(step))
    
    if num_halos <= 1.5:
        print("No halos in this snapshot.")
        return None
    
    # The path for the sim we want to analyze.
    LIST_PATH = '/scratch/hc2347/references/' + step.split('.')[-1] + '_filtered_halos.p'
    
#     if os.path.exists(LIST_PATH):
#         print('Found existing filtered halo list.')
#         halolist = pickle.load(open(LIST_PATH,'rb'))
#     else:
#         print('No halo list found. Creating new list.')
#         halolist = makeFilteredHalos(halo_catalogue)
#         pickle.dump(halolist, open(LIST_PATH,'wb'))

    halolist = makeFilteredHalos(halo_catalogue)
        pickle.dump(halolist, open(LIST_PATH,'wb'))
    
    # Distribute
    print('Spawning processes.')
    with concurrent.futures.ProcessPoolExecutor(max_workers=PROC_N) as executor:
        # result = executor.map(make_future_runner, range(1,15))
        result = executor.map(make_future_runner, halolist)
        flat_results = list(result)
    
    # Save it all
    if pickle_path is None:
            pickle_path = "/scratch/hc2347/pickles/60"
    if os.path.isdir(pickle_path) is False:
        os.mkdir(pickle_path)
    fname = os.path.join(pickle_path,'main_preload_z{:.3f}.p'.format(zred))
    print('Writing out {}'.format(fname))
    pickle.dump( flat_results, open( fname, 'wb' ))
    print('Done')

if __name__ == '__main__':
    
    '''
    Usage: python main/data/make_data.py /scratch/kld8/simulations/LRZ_Planck60 planck.new.hydro.60_600.00832
    '''
    
    path = sys.argv[1]
    step = sys.argv[2]
    
    # The path for the sim we want to analyze.
    SIM_FILE = os.path.join(path, step)
    LIST_PATH = '/scratch/hc2347/references/' + step.split('.')[-1] + '_filtered_halos.p'
    
    pynbody.config['threading'] = False
    
    sim = pynbody.load(SIM_FILE)
    halo_catalogue = sim.halos()

    zred = sim.properties['z']    
    num_halos = len(halo_catalogue)
    
    
    main(path, step)
