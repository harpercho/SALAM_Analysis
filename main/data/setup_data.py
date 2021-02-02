import os, glob

# Location of galaxies you want
sim_path = '/scratch/database/nihao/nihao_classic/'

# Folder that you want to put the galaxies in
parent_dir = '/scratch/hc2347/Planck60/'

# Decision whether to write out the list of galaxies
write_file = True

# Format of name of simulations, could be modified also
# e.g. 'g*12' to get all the 10^12 halos or g2.63e10 for a single halo
contained_files = glob.glob(os.path.join(sim_path,'*'))

# If writing file, give filename of the last folder, here is 'nihao_classic_dirs.txt'
if write_file:
    sim_type = sim_path.split('/')[-1]
    print('sim_type',sim_type,sim_path.split('/'))
    fname = '{}_dirs.txt'.format(sim_type)
    f = open(fname,'w')
    
# Loop through the matching halos
for file in contained_files:
    sim_dir = file
    halo_dir = halo_dir.split('/')[-1]
    # Following line is only to ignore files or directories, may be unecessary
    if '_' in halo_dir or '.p' in halo_dir or '.t' in halo_dir:
        continue
    print(halo_dir)
    
    if write_file:
        f.write('{}\n'.format(halo_dir))
​
    # Make parent directory and subdirectories for halos
    # Also, make your own subdirectories, here 'analysis', but maybe 'plots'
    # Or a list like for subfolder in ['analysis','plots','pickles']:
    try:
        os.makedirs(parent_dir)
    except OSError:
        pass
​
    source_path = os.path.join(parent_dir, halo_dir)
    try:
        os.makedirs(source_path)
    except OSError:
        pass
​
    source_path = os.path.join(parent_dir, halo_dir, 'analysis')
    try:
        os.makedirs(source_path)
    except OSError:
        pass
​
    sim_files = glob.glob(os.path.join(sim_dir,'*'))
    for sim_file in sim_files:
        source_path = sim_file
        sim_file = sim_file.split('/')[-1]
        out_path = os.path.join(parent_dir, halo_dir, sim_file)
        try:
            os.symlink(source_path, out_path)
        except OSError:
            pass
​
if write_file:
