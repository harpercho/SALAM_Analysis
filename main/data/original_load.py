import os
import numpy as np
import pickle
import pynbody
import pynbody.filt as filt
import multiprocessing
import concurrent.futures

def make_future_runner(idx):

    halo = halo_catalogue.load_copy(idx)
    #halo = halo_catalogue[idx]
    halo.physical_units()
    
    print("{}: Working on halo {}".format(multiprocessing.current_process().name, idx),flush=True)

    mvir=np.sum(halo['mass'].in_units('Msol'))
    mstar = np.sum(halo.star['mass'].in_units('Msol'))
    mgas = np.sum(halo.gas['mass'].in_units('Msol'))
    mdm = np.sum(halo.dm['mass'].in_units('Msol'))

    fifmyrf = filt.LowPass('age','10 Myr')
    sfr = np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / 1.e7

    fifmyrf = filt.LowPass('age','100 Myr')
    sfr_100 = np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / 1.e8
    
    print('Halo {} complete'.format(idx),flush=True)
    return {'mvir': mvir, 'mstar': mstar, 'mgas':mgas,'sfr':sfr,'sfr_100':sfr_100}


def main(path, step, pickle_path=None):
    
    print("Loading sim {}".format(simfile))
    
    if num_halos <= 1.5:
        print("No halos")
        return None
    
    LIST_PATH = '/scratch/hc2347/references/' + step.split('.')[-1] + '_filtered_halos.p'
    if os.path.exists(LIST_PATH):
        print('Found existing filtered halo list.')
        halolist = pickle.load(open(LIST_PATH,'rb'))
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        result = executor.map(make_future_runner, halolist)
        #result = executor.map(make_future_runner, range(1,3))
        flat_results = list(result)
        
    if pickle_path is None:
            pickle_path = "/scratch/hc2347/data/60/"
    if os.path.isdir(pickle_path) is False:
        os.mkdir(pickle_path)
    fname = os.path.join(pickle_path,'vol_halo_center_z{:.3f}.p'.format(zred))
    print('Writing out {}'.format(fname))
    pickle.dump( flat_results, open( fname, 'wb' ))
    print('Done')

if __name__ == '__main__':
    import sys
    path = sys.argv[1]
    step = sys.argv[2]
    pickle_path = None
    
    simfile = os.path.join(path, step)
    
    sim = pynbody.load(simfile)
    zred = sim.properties['z']
    halo_catalogue = sim.halos(dummy=True)
    num_halos = len(halo_catalogue)
    
    main(path, step, pickle_path)
