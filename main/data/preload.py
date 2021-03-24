import os
import sys
import numpy as np
import pickle
import pynbody
import pynbody.filt as filt
import multiprocessing
import concurrent.futures
from quantities import *
from haloprep import *
		
def make_future_runner(idx):
	# Function run by worker proceesses.
	
	print("{}: Working on halo {}/{}".format(multiprocessing.current_process().name, idx ,num_halos), flush = True)
	
	halo_dict = {}
	
	# Load a copy of the target halo.
	halo = halo_catalogue.load_copy(idx)

	# Center halo
	center(halo) #applies a transformation
	
	try:
		smallhalo = half_stellar_radius(halo)
	except Exception as e:
		print(e)
		print("Sticking with regular halo.")
		smallhalo = halo
	
	# Give physical units.
	halo.physical_units()
	
	coldgastemp = 1.5e4
	
	# From quantities.py. This file contains all the actual quantity obtaining utilities.
	# Each of these functions require a storage dictionary to update/add into.
	print("starting obtaining process for halo " + str(idx),flush=True)
	getParticleInfo(halo, halo_dict)
	getParticleInfo(smallhalo, halo_dict)
	getMasses(halo, halo_dict)
	getSFR(smallhalo, halo_dict, 10)
	getSFR(smallhalo, halo_dict, 100)
	getColdGas(halo, halo_dict, coldgastemp)
	getOxygenAbundance(smallhalo, halo_dict, coldgastemp)
	getStellarMetallicity(smallhalo, halo_dict, halo_dict['mstar'])
	
	# ADD YOUR OWN QUANTITIES HERE
	# For example:
	halo_dict['ID'] = idx
	halo_dict['mvir'] = np.sum(halo['mass'].in_units('Msol'))
	halo_dict['mstar_small'] = np.sum(smallhalo.s['mass'].in_units('Msol'))
	
	print('Halo ' + str(idx) + ' complete')
	
	return halo_dict

def main(path, step):
	
	PROC_N = 5 # number of worker processes
	
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
			

	# Distribute
	with concurrent.futures.ProcessPoolExecutor(max_workers=PROC_N) as executor:
		#result = executor.map(make_future_runner, range(15098,15105))
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
	Usage: python preload.py /scratch/kld8/simulations/LRZ_Planck60 planck.new.hydro.60_600.00832
	'''
	
	path = sys.argv[1]
	step = sys.argv[2]
	
	box = step.split('.')[3].split('_')[0]
	print(box)
	
	# The path for the sim we want to analyze.
	SIM_FILE = os.path.join(path, step)
	LIST_PATH = '/scratch/hc2347/references/' + box + '_' + step.split('.')[-1] + '_filtered_halos.p'


	pynbody.config['threading'] = False
	pynbody.config['number_of_threads'] = 1
	print("Main threads: " + str(pynbody.config['number_of_threads']), flush=True)
	
	sim = pynbody.load(SIM_FILE)
	try:
		halo_catalogue = sim.halos()
	except:
		print("No halo catalogue for " + step)
		sys.exit()
		
	halo_properties = sim.halos(dummy=True)

	zred = sim.properties['z']	  
	num_halos = len(halo_catalogue)
	
	# IMPORTANT! Change this output path to whatever you need
	OUTPUT = os.path.join("/scratch/hc2347/pickles/" + box,'oxhonly_{:.3f}.p'.format(zred))
	
	
	main(path, step)
