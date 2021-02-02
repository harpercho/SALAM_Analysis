import subprocess
import os, glob, pickle

def main(size=60):
    data_directory = "/scratch/kld8/simulations/LRZ_Planck{}".format(size)
    sims_with_halos = glob.glob(os.path.join(data_directory,'*.AHF_halos'))
    sims_with_halos.sort(reverse=True)

#     sims_with_halos = ['/scratch/kld8/simulations/LRZ_Planck60/planck.new.hydro.60_600.00448.z0.946.AHF_halos', 
#                        '/scratch/kld8/simulations/LRZ_Planck60/planck.new.hydro.60_600.00156.z2.994.AHF_halos']
    
    n_sims = len(sims_with_halos)

    for idx, sim in enumerate(sims_with_halos):
        iout = sim.split(".z")[0]
        file = iout.split("/")[-1] 
        print("Loading data for {}".format(file))

        gen_data = subprocess.run(["python","/scratch/hc2347/main/data/preload.py",data_directory,file])

#         stdout_path ="/scratch/hc2347/references/"+file+".txt"

#         os.makedirs(os.path.dirname(stdout_path), exist_ok=True)

#         with open(stdout_path,"wb") as f:
#             f.write(gen_data.stdout)
    
        print("Completed sim {}/{}".format(idx+1,n_sims))

if __name__ == '__main__':
    import sys
    size = sys.argv[1]
    main(size)
