import subprocess
import os, glob, pickle

# This is just a simple wrapper to run the preload script on all the available halo snapshots for SALAM.

def main(size=60):
    data_directory = "/scratch/kld8/simulations/LRZ_Planck{}".format(size)
    sims_with_halos = glob.glob(os.path.join(data_directory,'*.AHF_halos'))
    sims_with_halos.sort(reverse=True)
    sims_with_halos = sims_with_halos[1:]
#     sims_with_halos = ['/scratch/kld8/simulations/LRZ_Planck60/planck.new.hydro.60_600.00448.z0.946.AHF_halos', 
#                        '/scratch/kld8/simulations/LRZ_Planck60/planck.new.hydro.60_600.00156.z2.994.AHF_halos']
    print(sims_with_halos)
    n_sims = len(sims_with_halos)

    for idx, sim in enumerate(sims_with_halos):
        iout = sim.split(".z")[0]
        file = iout.split("/")[-1] 
        print("Loading data for {}".format(file))

        gen_data = subprocess.run(["python","/scratch/hc2347/main/data/preload.py",data_directory,file])

        stdout_path ="/scratch/hc2347/references/"+file+".txt"

        print("Completed sim {}/{}".format(idx+1,n_sims))

if __name__ == '__main__':
    import sys
    size = sys.argv[1]
    main(size)
