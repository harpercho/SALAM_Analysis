import pynbody
import numpy as np
import sys
from time import perf_counter

def create_coords(catalog):
	coords = np.zeros((len(catalog), 3))
	print(coords.shape)
	for i, halo in enumerate(catalog):
		coords[i][0], coords[i][1], coords[i][2] = halo.properties['Xc'], halo.properties['Yc'], halo.properties['Zc']
	
	return coords
		

def isIsolated(idx):
	halo = catalog[idx]
	x, y, z = halo.properties['Xc'], halo.properties['Yc'], halo.properties['Zc']
	dist_lim = halo.properties['Rvir']
	
	i = 1
	while i < idx:
		other_halo = coords[i]
		x_, y_, z_ = other_halo[0], other_halo[1], other_halo[2]
		dist = ((x_ - x)**2 + (y_ - y)**2 + (z_ - z)**2)**0.5
		
		if dist < dist_lim:
			return False
		i += 1
	return True

if __name__ == '__main__':

	#path = sys.argv[1]
	path = "/scratch/kld8/simulations/LRZ_Planck60/planck.new.hydro.60_600.00832"
	
	s = pynbody.load(path)
	catalog = s.halos()
	
	coords = create_coords(catalog)
	print("The time taken to preload was ", end - start)
	
	idxs = np.random.randint(1, len(catalog), 20)
	print("Is halo {} isolated?".format(idx))
	start = perf_counter()
	for idx in idxs:
		print(isIsolated(idx))
	end = perf_counter()

	print("The time taken for the loop is {} s".format(end-start))
