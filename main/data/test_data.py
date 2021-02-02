import os
import numpy as np
import pickle
import pynbody
import pynbody.filt as filt
import multiprocessing as mp
import concurrent.futures

def f(idx):
    print(multiprocessing.current_process())
    halo = halo_catalogue[idx]
    print(halo)

def main(path, step):
    print("Started main")
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        result = executor.map(f, range(1,15))
        flat_results = list(result)
    print(flat_results)
    

if __name__ == '__main__':
    import sys
    path = sys.argv[1]
    step = sys.argv[2]
    
    # The path for the sim we want to analyze.
    SIM_FILE = os.path.join(path, step)
    LIST_PATH = '/scratch/hc2347/references/' + step.split('.')[-1] + '_filtered_halos.p'
    
    sim = pynbody.load(SIM_FILE)
    halo_catalogue = sim.halos()

    zred = sim.properties['z']    
    num_halos = len(halo_catalogue)
    
    
    main(path, step)
