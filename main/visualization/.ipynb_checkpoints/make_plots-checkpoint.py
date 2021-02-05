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
    
def plot_SMF():
    ax.set_yscale("log")
    plot_mf(ax,filter_list(entry["mstar"],10**7, 10**12),50,entry["zred"])
    
    plt.legend(frameon=False)
    ax.legend()
    
    plt.savefig('/scratch/hc2347/reports/60/center_smf60.png')
    
def plot_Moster():
    from pynbody.plot.stars import moster
    xlabel = '$\\rm{log_{10}}(M_{200}/\\rm M_{\odot})$'
    ylabel = '$\\rm{log_{10}}(M_{\star}/\\rm M_{\odot})$'
    label = '$z \;=\; {}$'.format(zred)
    c = 'steelblue'
    
    #Moster
    xmasses = np.logspace(np.log10(min(entry['mvir'])),1+np.log10(max(entry['mvir'])),20)
    ystarmasses, errors = moster(xmasses,float(entry["zred"]))
    ax.plot(np.log10(xmasses),np.log10(np.array(ystarmasses)),color="Black",label="Moster 2013")
    ax.plot(np.log10(xmasses),np.log10(np.array(ystarmasses)/np.array(errors)), linestyle='dashed', color = 'grey')
    ax.plot(np.log10(xmasses),np.log10(np.array(ystarmasses)*np.array(errors)), linestyle='dashed', color = 'grey')
    ax.fill_between(np.log10(xmasses),np.log10(np.array(ystarmasses)/np.array(errors)),y2=np.log10(np.array(ystarmasses)*np.array(errors)),
                    color='grey', alpha=0.2)
    

