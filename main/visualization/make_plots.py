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
    """Plot the stellar mass function and save reports directory."""
    
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
        
    #Axes settings
    ax.set_yscale('log')
    ax.legend(frameon=False)
    ax.tick_params(direction='in', which='both')
    plt.tight_layout()
    
    plt.savefig("/scratch/hc2347/reports/60/SMF/z{}.png".format(z))

    
def plot_coldgas(mgas,mstar):
    
    mstar, mgas = do_filter(mstar, mgas)
    
    fg = mgas/mstar
    
    x = np.log10(mstar)
    y = np.log10(fg)

    
    fig, ax = plt.subplots(figsize = (9,7))
    
    print(str(len(mstar)) + ' number of points.')
    
    hb = ax.hexbin(x,y,gridsize= 50, bins='log',cmap='Blues')
    cb = fig.colorbar(hb, ax=ax)
    
    peebles = np.genfromtxt('/scratch/hc2347/references/obs/Peebles_2014_coldgas.csv',unpack=True,skip_header=2,delimiter=',')
    logmstar_pb = peebles[0]
    median_pb = peebles[1]
    sixteen_pb = peebles[2]
    eightyfour_pb = peebles[3]
    
    nihao_x = nihao('Mstar',0)
    nihao_y = nihao('Mcool',0)/nihao('Mstar',0)
    print(nihao_x[:15])
    print(nihao_y[:15])
        
    ax.scatter(np.log10(nihao_x), np.log10(nihao_y), color='maroon', label = "NIHAO Classic")
    
    ax.plot(logmstar_pb, np.log10(median_pb), c='black', label = "Peeples 2014")
    ax.plot(logmstar_pb, np.log10(sixteen_pb), linestyle='dashed',color='grey')
    ax.plot(logmstar_pb, np.log10(eightyfour_pb), linestyle='dashed',color='grey')
    ax.fill_between(logmstar_pb, np.log10(sixteen_pb), y2 = np.log10(eightyfour_pb), alpha = 0.1, color='grey')
    ax.set_ylabel('$log_{10}F_g \equiv log_{10}(M_{gas}/M_{*})$',fontsize=12)
    ax.set_xlabel('$ log M_{*}/M_\odot$',fontsize=12)
    
    ax.set_xlim(7,13)
    ax.set_ylim(-2,2)
    #ax.set_yscale('log')
    ax.legend()

    
    plt.title("Cold Gas Fraction")
    plt.savefig("/scratch/hc2347/reports/cold_gas.png")

    
def plot_oxh(oxh,mstar):
    
    tremonti = np.genfromtxt('/scratch/hc2347/references/obs/Tremonti_2004_mzr.csv',unpack=True,skip_header=2,delimiter=',')
    
    median_tr = tremonti[3]
    sixteen_tr = tremonti[2]
    eightyfour_tr = tremonti[4]
    
    logmstar_tr =tremonti[0]
    
    fig, ax = plt.subplots(figsize = (9,7))
    
    x = np.log10(mstar)
    y = oxh
    
    x,y = do_filter(x,y)
    
    hb = plt.hexbin(x, y, gridsize= 100, bins='log',cmap='Blues')
    plt.colorbar(hb)
    
    ax.plot(logmstar_tr, median_tr, c='black', label = "Tremonti 2004 fit")
    ax.plot(logmstar_tr, sixteen_tr, linestyle='dashed',color='grey')
    ax.plot(logmstar_tr, eightyfour_tr, linestyle='dashed',color='grey')
    ax.fill_between(logmstar_tr, sixteen_tr, y2 = eightyfour_tr, alpha = 0.1, color='grey')
    ax.set_ylabel('$12+\log_{10}(O/H))$',fontsize=12)
    ax.set_xlabel('$ log M_{*}/M_\odot$',fontsize=12)
    
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()

    ax.axis([xmin, xmax, ymin, ymax])
    
    nihao_x = nihao('Mstar',0)
    nihao_y = nihao('oxh',0)
        
    ax.scatter(np.log10(nihao_x), nihao_y, color='maroon', label = "NIHAO Classic")
    
    ax.legend()
    ax.set_ylim(6,10)
    #ax.set_xlim(8.5,12)
    
    plt.savefig("/scratch/hc2347/reports/60/CenterCold_MZR_Oxh.png")
