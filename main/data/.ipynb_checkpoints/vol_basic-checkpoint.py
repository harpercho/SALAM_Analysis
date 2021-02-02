import os, glob
import numpy as np
import pickle
import pynbody
import pynbody.filt as filt

#pynbody.config['number_of_threads'] = 1 - set depending on number of processes

def halo_ix_init(h, shuffle=False):
    '''
    Return list of indicies to halos, which can be scattered with MPI
    shuffle (bool) - whether to shuffle the array
    '''
    halo_ix = list(range(1,len(h)+1))
    if shuffle:
        import random
        random.shuffle(halo_ix)
    return halo_ix

def main(path, iout, pickle_path=None):
    from parallel import mpi 

    sim_file = os.path.join(path,iout)
    mpi.msg('Loading sim_file {}'.format(sim_file))
    s = pynbody.load(sim_file)
    #s.set_nproc(1)
    h = s.halos()
    if len(h) < 1.5:
        mpi.msg('No halos...')
        return None
    zred = s.properties['z']
    
    halo_ix = None
    if mpi.host:
        halo_ix = halo_ix_init(h, shuffle=True)

    dest = {}
    count = 0
    for i, sto in mpi.piter(halo_ix, storage=dest):
        halo = h.load_copy(i)
        halo.physical_units()
        try:
        #    pynbody.analysis.angmom.faceon(halo)
            pynbody.analysis.halo.center(halo,mode='pot')
        except:
            count += 1
            mpi.msg('Cannot center halo'.format(count))

        mpi.msg('Working on halo {}'.format(i))

        mstar = np.sum(halo.star['mass'].in_units('Msol'))
        mgas = np.sum(halo.gas['mass'].in_units('Msol'))

        fifmyrf = filt.LowPass('age','10 Myr')
        sfr = np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / 1.e7
        
        fifmyrf = filt.LowPass('age','100 Myr')
        sfr_200 = np.sum(halo.star[fifmyrf]['mass'].in_units('Msol')) / 1.e8
        
        '''try: - preferred method at late times, and I think you have the work around
            init_mstar = s.star['massform'].in_units('Msol')[0]
            fifmyrf = filt.LowPass('age','30 Myr')
            sfr_i = len(halo.star[fifmyrf])*init_mstar / 3.e7
        except:
            init_mstar = False
            sfr_i = 0.
        mpi.msg('Other sfr {}'.format(sfr_i))'''
        mvir=np.sum(halo['mass'].in_units('Msol'))
        
        sto.idx = len(halo)
        sto.result = {'Mstar' : mstar, 'SFR' : sfr, 'SFR_200' : sfr_200, \
                          'Mvir' : mvir, 'Mgas' : mgas}
    if mpi.host:
        if pickle_path is None:
            pickle_path = os.path.join(path,'analysis')
        if os.path.isdir(pickle_path) is False:
            os.mkdir(pickle_path)
        fname = 'vol_halo_av_z{:.3f}.p'.format(zred)
        mpi.msg('Writing out {}'.format(fname))
        pickle.dump( mpi.unpack(dest), open( fname, 'wb' ) )
        mpi.msg('Done')


if __name__ == "__main__":
    import sys
    path = sys.argv[1]
    iout = sys.argv[2]
    pickle_path = None
    #if len(sys.argv) > 4:
     #   pickle_path = sys.argv[4]

    main(path, iout, pickle_path)
