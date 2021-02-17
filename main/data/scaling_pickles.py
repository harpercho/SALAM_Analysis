# Simple script to dump a bunch of halo calculations into pickles
# Each halo gets a pickle for the chosen snaps, then all info dumped
# The annoying try/excepts are for when there are no star particles or halos.

import os, glob, pickle
import numpy as np
import pynbody
import pynbody.filt as filt

#sim_names = glob.glob(os.path.join('..','g*'))
#sim_names = ['/scratch/kld8/simulations/g9.59e10']
#bad_halos = ['g1.05e13','g8.94e12']#['g4.55e08','g4.00e10','g1.07e09','g2.21e09','g1.15e10']
sim_file = open('nihao_classic_dirs.txt','r')
sim_names = sim_file.read().split('\n')[:-2]
sim_file.close()
print(sim_names)

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
    halo_dict_file = os.path.join('/scratch/hc2347/pickles',sim_name,'analysis','pickle_Harper.p')
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
        halo_dict = {}
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
                    
            mstar = np.sum(halo.s['mass'].in_units('Msol'))
            mdm = np.sum(halo.d['mass'].in_units('Msol'))

            mgas = np.sum(halo.g['mass'].in_units('Msol'))
            temp = 1.5e4
            coolgasf = filt.And(filt.LowPass('temp',temp),filt.HighPass('rho','0.03 m_p cm^-3'))
            mcool = np.sum(halo.g[coolgasf]['mass'].in_units('Msol'))
            #print('M_star {} and M_gas {}'.format(mstar, mgas))

            fifmyrf = filt.LowPass('age','10 Myr')
            sfr_10 = np.sum(halo.s[fifmyrf]['mass'].in_units('Msol')) / 1.e7

            fifmyrf = filt.LowPass('age','100 Myr')
            sfr_100 = np.sum(halo.s[fifmyrf]['mass'].in_units('Msol')) / 1.e8

            mvir=np.sum(halo['mass'].in_units('Msol'))
            
            z_star = np.sum(halo.s['mass']*halo.s['metals'])/mstar

            oxh = np.log10(np.sum(halo.g['OxMassFrac'])/ \
                               (16.*np.sum(halo.g['hydrogen']))) + 12.
            
            halo_dict[zred] = {}
            halo_dict[zred]['Mstar'] = mstar
            halo_dict[zred]['Mdm'] = mdm
            halo_dict[zred]['Mgas'] = mgas
            halo_dict[zred]['SFR_10'] = sfr_10
            halo_dict[zred]['SFR_100'] = sfr_100
            halo_dict[zred]['Mvir'] = mvir
            halo_dict[zred]['Mcool'] = mcool
            halo_dict[zred]['z_star'] = z_star
            halo_dict[zred]['oxh'] = oxh
        gal_dict[halo_name] = halo_dict
        pickle.dump( halo_dict, open( halo_dict_file, 'wb') )
    #except:
     #   print('Bad_halo {}'.format(halo_name))
     #   pass
pickle.dump( gal_dict, open('pickle_NIHAO.p','wb') )

