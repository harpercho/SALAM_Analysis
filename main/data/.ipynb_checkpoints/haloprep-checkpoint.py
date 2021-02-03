def isHostHalo(h):
    return h.properties['hostHalo'] <= 0
    
def isBigHalo(h,cut):
    return len(h) > cut
    
def makeFilteredHalos(h):
    
    '''
    Creates a index list of halos that pass our filter.
    Default is 100 particles (generic) and must be hosthalo.
    
    If WriteFile is true, it writes this into a pickle file in ../references.
    '''
    
    # Set minimum number of particles
    P_MIN = 100
    
    halolist = []
    failed = 0
    
    for i,halo in enumerate(h):
        if isHostHalo(halo) and isBigHalo(halo,P_MIN):
            halolist.append(i+1)
            print("Appended " + str(i))
        else:
            print("Removed " + str(i))
            failed +=1
            
    print('Removed ' + str(failed) + 'halos')
    
    return halolist

def center(halo):
    try:
        # Make vel=True only if velocity is useful in the quantities being obtained.
        # By default uses 'ssc' mode according to pynbody.config
        centering = pynbody.analysis.halo.center(halo, vel=False)
    except Exception as e:
        print("Could not center: {}".format(e), flush=True)
    
def radius_cut(halo, idx):
    try:
        RVIR = pynbody.array.SimArray(np.max(halo['r'].in_units('kpc')),'kpc')
    except Exception as e:
        print("No RVIR. Trying AHF method.")
        print(e)
        try:
            RVIR = halo_properties[idx].properties['Rvir']
        except:
                print("No RVIR.")
                return halo
                
    diskf = filt.Sphere(str(RVIR*0.2) +' kpc') #20% Virial Radius
    
    return halo[diskf]