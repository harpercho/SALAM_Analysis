import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import matplotlib.pyplot as plt
import sys, os, glob, pickle, struct
import plot_tools

def plot_HMF(simpath, entry, ax=None):
    """Plot the halo mass function. If there are no axes given, make a figure and axes"""
    
    z = entry["zred"]
    
    if ax == None:
        fig, ax = plt.subplots(figsize=(7,5), dpi = 300)
        
    s = pynbody.load(simpath)
    
    x, hmf = plot_tools.plot_mf(ax,plot_tools.filter_list(entry["mvir"],10**10,10**15),50,z)
    ax.plot(x, hmf, label = "SALAM z{}".format(entry['zred']))
    stms, stsig, stmf = pynbody.analysis.halo_mass_function(s)
    ax.semilogy(np.log10(stms),stmf,label="Sheth-Tormen", color = "#808080", linestyle="--")
    
    ax.set_xlabel('$\\rm{log_{10}}(M_{vir}/\\rm M_{\odot})$')
    ax.set_ylabel('number density [Mpc$^{-3}$]')
    
    ax.legend(frameon=False)
    ax.tick_params(direction='in', which='both')
    plt.tight_layout()
    
    plt.savefig("/scratch/hc2347/reports/60/HMF/z{}.png".format(z))
    
def plot_SMF(entry, bins):
    """Plot the stellar mass function and save reports directory."""
    
    z = entry["zred"]
    stellar_mass_density = plot_tools.filter_list(entry["mstar"],10**7, 10**12)
        
    fig, ax = plt.subplots(figsize=(7,5), dpi = 150)
    
    # Plot Stellar Mass Function from density
    x, smf = plot_tools.plot_mf(ax, stellar_mass_density, bins, entry["zred"])
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
    
    plt.savefig("/scratch/hc2347/reports/60/SMF/z{}_bins{}.png".format(z, bins))

    
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
    
    fig, ax = plt.subplots(figsize = (9,7))

    # Add tremonti observations
    tremonti = np.genfromtxt('/scratch/hc2347/main/visualization/obs/Tremonti_2004_MZR_z0.csv',unpack=True,skip_header=2,delimiter=',')
    median_tr = tremonti[3]
    sixteen_tr = tremonti[2]
    eightyfour_tr = tremonti[4]
    logmstar_tr =tremonti[0]

    x = np.log10(mstar)
    y = oxh
    print(len(y))
    x,y = plot_tools.do_filter(x,y)
    
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
    
    nihao_x = plot_tools.nihao('mstar',0)
    nihao_y = plot_tools.nihao('oxh',0)
        
    ax.scatter(np.log10(nihao_x), nihao_y, color='maroon', label = "NIHAO Classic")
    
    ax.legend()
    ax.set_ylim(6,10)
    ax.set_xlim(7,12)
    plt.show()
    plt.savefig("/scratch/hc2347/reports/60/CenterCold_MZR_Oxh.png")
    
def plot_Moster(entry):
    from pynbody.plot.stars import moster

    z = entry["zred"]
    xlabel = '$\\rm{log_{10}}(M_{200}/\\rm M_{\odot})$'
    ylabel = '$\\rm{log_{10}}(M_{\star}/\\rm M_{\odot})$'
    label = '$z \;=\; {}$'.format(z)
    
    fig, ax = plt.subplots(figsize=(7,5), dpi = 300)
    
    x = np.log10(np.array(entry["mvir"]))
    y = np.log10(np.array(entry["mstar"]))
    x, y = plot_tools.do_filter(x,y)
    ax.hexbin(x, y, gridsize=50,cmap='Blues', bins='log')
    
    # Moster
    xmasses = np.logspace(np.log10(min(entry['mvir'])),1+np.log10(max(entry['mvir'])),20)
    ystarmasses, errors = moster(xmasses,float(entry["zred"]))
    ax.plot(np.log10(xmasses),np.log10(np.array(ystarmasses)),color="Black",label="Moster 2013")
    ax.plot(np.log10(xmasses),np.log10(np.array(ystarmasses)/np.array(errors)), linestyle='dashed', color = 'grey')
    ax.plot(np.log10(xmasses),np.log10(np.array(ystarmasses)*np.array(errors)), linestyle='dashed', color = 'grey')
    ax.fill_between(np.log10(xmasses),np.log10(np.array(ystarmasses)/np.array(errors)),y2=np.log10(np.array(ystarmasses)*np.array(errors)),
                    color='grey', alpha=0.2)
    
    # NIHAO
    nihao_mstar = plot_tools.nihao('Mstar',0)
    nihao_mhalo = plot_tools.nihao('Mvir',0)
    ax.scatter(np.log10(nihao_mhalo), np.log10(nihao_mstar), color='maroon', label = 'Nihao Classic')

    
    # Axes params
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    ax.legend(frameon=False)
    ax.tick_params(direction='in', which='both')
    ax.set(xlabel=xlabel, ylabel=ylabel)
    plt.tight_layout()
    
    plt.savefig("/scratch/hc2347/reports/60/MstarMhalo/Moster_z0.png")

def plot_SFR(entry):
    def salim(log_mstar):
        sfr = 5.96e-11*10**((-1.35 + 1)*(log_mstar - 11.03))*np.exp(-10**(log_mstar - 11.03))
        return sfr
    whitacker = np.genfromtxt('/scratch/hc2347/main/visualization/obs/Whitacker_2014_SMF_z0.5.csv',unpack=True,skip_header=1,delimiter=',')

    xlabel = '$\\rm{log_{10}}(M_{\star}/\\rm M_{\odot})$'
    ylabel = '$\\rm{log_{10}(SFR/M_{\odot}\,yr^{-1})}$'
    
    fig, ax = plt.subplots(figsize=(7,5), dpi = 300)
    
    xmin = 7
    xmax = 12
    ymin = -1
    ymax = 3
    
    x = np.log10(np.array(entry["mstar"]))
    y = np.log10(np.array(entry["sfr_100"]))
    #x, y = plot_tools.do_filter(x,y)
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))

    ax.hexbin(x, y,  gridsize = 100, bins = 'log', cmap = 'Blues')

    nihao_x = plot_tools.nihao('Mstar',0)
    nihao_y = plot_tools.nihao('SFR_100',0)
    ax.scatter(np.log10(nihao_x), np.log10(nihao_y), color='maroon', label="NIHAO Classic")
    #ax.plot(x,np.log10(salim(x)))
    ax.plot(whitacker[0], whitacker[1], marker='s', linestyle='none', color = 'orange', label='Whitacker 2014')
    ax.errorbar(whitacker[0], whitacker[1], yerr = whitacker[2], fmt='none', color = 'orange')
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(direction='in', which='both')
    ax.legend(frameon=False)
    plt.tight_layout()
    
    plt.savefig('/scratch/hc2347/reports/60/sfr_100_v1.png')
    
def plot_stellar_metallicity(z_star, m_star):

    z_sol = 0.013 # primordial Solar metallicity


    x = np.log10(m_star)
    y = np.log10(z_star/z_sol)
    gallazzi = np.genfromtxt('/scratch/hc2347/main/visualization/obs/Gallazzi_2005_SMZ_z0.csv',unpack=True,skip_header=2,delimiter=',')
    logmstar_tr = gallazzi[0]

    median_tr = gallazzi[1]
    sixteen_tr = gallazzi[2]
    eightyfour_tr = gallazzi[3]

    fig, ax = plt.subplots(figsize=(8,6))
    
    plt.hexbin(x,y,gridsize=100,bins='log',cmap="Blues")
#     cb = plt.colorbar()
#     cb.set_label("counts")

    nihao_x = np.array(plot_tools.nihao('mstar',0))
    nihao_y = np.array([np.sum(z_vals) for z_vals in plot_tools.nihao('z_star',0)])
        
    print(nihao_x[:10])
    print(nihao_y[:10])
        
    ax.scatter(np.log10(nihao_x), np.log10(nihao_y/z_sol), color='maroon', label = "NIHAO Classic")
    
    
#     ax.plot(logmstar_tr, median_tr, c='black',label = "Gallazzi 2005")
#     ax.plot(logmstar_tr, sixteen_tr, linestyle='dashed',color='grey')
#     ax.plot(logmstar_tr, eightyfour_tr, linestyle='dashed',color='grey')
    
#     ax.fill_between(logmstar_tr, sixteen_tr, y2 = eightyfour_tr, alpha = 0.2, color='grey' )
#     ax.set_ylabel('$log(Z_{*}/Z_\odot)$',fontsize=12)
#     ax.set_xlabel('$ log M_{*}/M_\odot$',fontsize=12)
    ax.legend()

    ax.set_ylim(-2,0.5)
    ax.set_xlim(7,12)
    
    plt.savefig("/scratch/hc2347/reports/60/Center_Gallazzi_SMZR.png")