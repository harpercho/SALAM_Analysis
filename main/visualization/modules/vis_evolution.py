import numpy as np
import pynbody
from modules.plot_tools import *
from modules.add_obs import *
import matplotlib.pyplot as plt
import sys, os, glob, pickle, pylab as plt, struct

DATA_PATH = "/scratch/hc2347/data/60/vol_halo_z" # The assumption being that all relevant volume pickles are of this format
SAVE_PATH = "/scratch/hc2347/reports/evolution"

# LOAD relevant pickles

pickle_files = glob.glob(DATA_PATH+"*.p")
pickle_files.sort()

all_dicts = [] # List to store dictionaries for each zred
zred = [] # List to store the redshifts

for pickle in pickle_files:
    data = load_halos_pickle(pickle)
    z = file.split("z")[1][:-2]
    
    all_dicts.append(data)
    zred.append(zred)

# --- Stellar Mass Function ---

fig, axs = plt.subplots(2, 4, sharey=True,sharex=True, figsize=(15,8))
fig.subplots_adjust(hspace=0,wspace=0)
fig.suptitle("SMF",fontsize=12,fontweight="bold", y=0.90)

ax_last = axs[1][3]


fig.text(0.5, 0.08, 'log$_{10}$(M)', ha='center',fontsize=15)
fig.text(0.08, 0.5, 'dN / dlog$_{10}$(M)', va='center', rotation='vertical',fontsize=15)

z_simple = [0.2,1,3,4,5,6,10,10]
import add_obs

for idx, entry in enumerate(entries):
    
    PER_ROW = 4
    y = int(idx/PER_ROW)
    x = idx % PER_ROW

    ax = axs[y][x]

    ax.set_yscale("log")
    
    plot_mf(ax,filter_list(entry["mstar"]),50,entry["zred"])
    
    add_obs.SMF_Moustakas(ax,z_simple[idx],color='#696969')
    add_obs.SMF_Song(ax,z_simple[idx])
    add_obs.SMF_Duncan(ax,z_simple[idx], color = "#BBBBBB")
    add_obs.SMF_Mortlock(ax, z_simple[idx],continuous = True)
    
    ax.legend()
    
    plot_mf(ax_last,filter_list(entry["mstar"]),50,entry["zred"])
    ax_last.legend()

plt.savefig("/scratch/hc2347/reports/SMF_60.png")


