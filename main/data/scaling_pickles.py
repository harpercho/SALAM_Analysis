# Simple script to dump a bunch of halo calculations into pickles
# Each halo gets a pickle for the chosen snaps, then all info dumped
# The annoying try/excepts are for when there are no star particles or halos.

import os, glob, pickle
import numpy as np
import pynbody
import pynbody.filt as filt
import haloprep
from quantities import *

#sim_names = glob.glob(os.path.join('..','g*'))
#sim_names = ['/scratch/kld8/simulations/g9.59e10']
#bad_halos = ['g1.05e13','g8.94e12']#['g4.55e08','g4.00e10','g1.07e09','g2.21e09','g1.15e10']
sim_file = open('/scratch/hc2347/main/data/nihao_classic_dirs.txt','r')
sim_names = sim_file.read().split('\n')[:-2]
sim_file.close()
print(sim_names)

def isPolluted(halo):
    dm_part = halo.d['mass'].in_units('Msol')
    M_res = dm_part.min()
    Npoll = len(dm_part[dm_part > M_res])
    return Npoll > 0.10*len(dm_part)

try:
    gal_dict = pickle.load(open('pickle_NIHAO.p','rb'))
    existing_keys = gal_dict.keys()

except:
    gal_dict = {}
    existing_keys = [] #bad_halos

for sim_name in sim_names:
    halo_name = sim_name #.split('/')[-1]
    if halo_name in existing_keys:
        print('Skipping halo {}'.format(halo_name))
        continue
    halo_dict_file = os.path.join('/scratch/hc2347/pickles/nihao/',sim_name,'pickle_run0.p')
    try:
        halo_dict = pickle.load(open(halo_dict_file,'rb'))
        gal_dict[halo_name] = halo_dict
        print('Halo pickle exists {}'.format(halo_dict_file))
    except:
        # These next two lines are because I wanted to specific redshifts
        sim_files = glob.glob(os.path.join('/scratch/kld8/simulations',sim_name,'g*00112*.HeII'))
        sim_files += glob.glob(os.path.join('/scratch/kld8/simulations',sim_name,'g*01024*.HeII'))
        sim_files = [x[:-5] for x in sim_files]
        print('List of outputs to analyze: {}'.format(sim_files))

        print('Doing halo {}'.format(halo_name))
        halo_dict_big = {}
        for sim_file in sim_files:
            print('Sim file name {}'.format(sim_file))
            s = pynbody.load(sim_file)
            zred = s.properties['z']
        
            h = s.halos()
            try:
                init_mstar = s.star['massform'].in_units('Msol')[0]
            except:
                init_mstar = False
            if len(h) < 1.5:
                continue
            else:
                try: 
                    pynbody.analysis.angmom.faceon(h[1])
                    s.physical_units()
                    halo = h[1]
                except:
                    print('Failed on halo {}'.format(halo_name))
                    s.physical_units()
                    halo = h[1]    
            curr_idx = 1
            while isPolluted(halo):
                curr_idx += 1
                halo = h[curr_idx]
                
            haloprep.center(halo)
            try:
                smallhalo = haloprep.half_stellar_radius(halo)
            except Exception as e:
                print(e)
                RVIR = pynbody.array.SimArray(np.max(halo['r'].in_units('kpc')),'kpc')
                diskf = filt.Sphere(str(RVIR*0.2) +' kpc')
                smallhalo = halo[diskf]
            
            halo_dict_big[zred] = {}            
            halo_dict = halo_dict_big[zred]
            
            coldgastemp = 1.5e4
            getParticleInfo(halo, halo_dict)
            getMasses(halo, halo_dict)
            getSFR(smallhalo, halo_dict, 10)
            getSFR(smallhalo, halo_dict, 100)
            getColdGas(halo, halo_dict, coldgastemp)
            getOxygenAbundance(halo, halo_dict, coldgastemp)
            getStellarMetallicity(smallhalo, halo_dict, halo_dict['mstar'])
            

        gal_dict[halo_name] = halo_dict_big
        
        if not os.path.exists(halo_dict_file):
            os.makedirs(halo_dict_file.split('pickle_')[0])
        pickle.dump( halo_dict, open(halo_dict_file, 'wb') )
    #except:
     #   print('Bad_halo {}'.format(halo_name))
     #   pass
pickle.dump( gal_dict, open('/scratch/hc2347/pickles/nihao/pickle_NIHAO.p','wb') )

