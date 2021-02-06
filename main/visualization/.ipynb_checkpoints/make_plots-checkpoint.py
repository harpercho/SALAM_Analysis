def plot_HMF(simpath, entry, ax=None):
    """Plot the halo mass function. If there are no axes given, make a figure and axes"""
    if ax = None:
        fig, ax = plt.subplots()
        
    s = pynbody.load(simpath)
    
    plot_mf(ax,filter_list(entry["mvir"],10**10,10**15),100,entry["zred"])
    stms, stsig, stmf = pynbody.analysis.halo_mass_function(s)
    ax.semilogy(np.log10(stms),stmf,label="Sheth-Tormen", color = "#808080")
    ax.legend()
    
    plt.savefig('/scractch/hc2347/rep')
    
def plot_SMF(entry):
    z = entry["zred"]
    stellar_mass_density = plot_tools.filter_list(entry["mstar"],10**7, 10**12)
        
    fig, ax = plt.subplots(figsize=(7,5), dpi = 150)
    
    # Plot Stellar Mass Function from density
    x, smf = plot_tools.plot_mf(ax, stellar_mass_density, 50, entry["zred"])
    ax.plot(x, smf, label = 'SALAM z{}'.format(z))
    ax.set_ylabel('number density [Mpc $^{-3}$]')
    ax.set_xlabel('log$_{10}$(M) [M$_\odot$]')

    
    # Add observational results
    if z < 1:
        baldry = np.genfromtxt('obs/Baldry_2012_SMF_z0.csv',unpack=True,skip_header=1,delimiter=',')
        moustakas = np.genfromtxt('obs/Moustakas_2013_SMF_z0.csv', unpack=True, skip_header=1, delimiter=',')
        ax.plot(baldry[0], baldry[1]*10e-4, marker = '.', label='Baldry 2012', linestyle='None')
        ax.plot(moustakas[0], 10**moustakas[1], marker = '+', label='Moustakas 2013', linestyle='None')
        ax.errorbar(baldry[0], baldry[1], fmt='none')

    ax.set_yscale('log')
    ax.legend(frameon=False)
    ax.tick_params(direction='in', which='both')
    plt.tight_layout()
    
    plt.savefig("/scratch/hc2347/reports/60/SMF/z{}.png".format(z))