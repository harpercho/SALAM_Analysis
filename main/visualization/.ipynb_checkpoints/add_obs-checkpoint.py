import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
import matplotlib.pyplot as plt
import sys, os, glob, pickle

def SMF_Mortlock(ax, z, continuous=False):
    ''' 
    Mortlock 2013
    Only contains data points for z 1 2 and 3
    '''
    if z not in [1,2,3]:
        return None

    if z == 1:
        mstar = 10**np.array([8.6,8.9,9.1,9.4,9.6,9.9,10.1,10.4,10.6,10.9,11.1,11.4,11.6,11.9])
        y = np.array([-1.79,-1.83,-1.99,-2.06,-2.31,-2.40,-2.48,-2.60,-2.67,-2.55,-2.89,-3.17,-3.47,-4.17])
        yerr = np.array([0.12,0.12,0.11,0.12,0.13,0.12,0.13,0.13,0.13,0.14,0.15,0.15,0.18,0.24])
    
    if z == 2:
        mstar = 10**np.array([9.4,9.6,9.9,10.1,10.4,10.6,10.9,11.1,11.4,11.6,11.9])
        y = np.array([-2.17,-2.15,-2.21,-2.32,-2.62,-2.70,-2.89,-3.45,-3.64 ,-3.98,-4.75])
        yerr = np.array([0.12,0.12,0.12,0.13,0.13,0.13,0.13,0.42,0.4,0.54,1.81])
    
    if z == 3:
        mstar = 10**np.array([9.6,9.9,10.1,10.4,10.6,10.9,11.1,11.4,11.6,11.9,12.1])
        y = np.array([-2.45,-2.56,-2.62,-2.88,-2.98,-3.14,-3.96,-3.78,-3.96,-3.96,-3.96])
        yerr = np.array([0.12,0.13,0.13,0.14,0.14,0.16,0.20,0.17,0.19,0.19,0.22])

    if continuous:
        y_upper = 10**np.array(y + yerr)
        y_lower = 10**np.array(y - yerr)
        ax.plot(np.log10(mstar), y, color="black")     
        ax.fill_between(np.log10(mstar), y_upper,  y_lower,
                         label='Mortlock 2013', color='grey', alpha=0.5)
        return None
    
    ax.plot(np.log10(mstar), 10**y, 'o')
    return None

def SMF_Song(ax, z, color='steelblue'):
    print("started")
    '''
    Song 2016
    High Redshift values only  4 < z < 8
    '''
    if z not in [4,5,6,7,8]:
        return None
    
    if z == 4:
        mstar = 10**np.array([7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25,10.75, 11.25])
        y = np.array([-1.57,-1.77,-2.00,-2.22,-2.52,-2.91,-3.37,-4.00,-4.54])
        y_upper = 10**np.array(y + np.array([0.21,0.15,0.13,0.09,0.09,0.12,0.09,0.20,0.34]))
        y_lower = 10**np.array(y - np.array([0.16,0.13,0.10,0.09,0.09,0.05,0.12,0.25,0.55]))

    if z == 5:
        mstar = 10**np.array([7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25,10.75, 11.25])
        y = np.array([-1.47, -1.72, -2.01, -2.33, -2.68, -3.12, -3.47, -4.12, -4.88])
        y_upper = 10**np.array(y + np.array([0.24,0.20,0.16,0.15,0.07,0.09,0.16,0.25,0.40]))
        y_lower = 10**np.array(y - np.array([0.21,0.20,0.16,0.10,0.14,0.11,0.14,0.38,0.61]))

    if z == 6:
        mstar = 10**np.array([7.25, 7.75, 8.25, 8.75, 9.25, 9.75,10.25])
        y = np.array([-1.47,-1.81,-2.26,-2.65,-3.14,-3.69,-4.27])
        y_upper = 10**np.array([-1.12, -1.58, -2.05, -2.50, -3.02, -3.57,-3.87])
        y_lower = 10**np.array([-1.79, -2.09, -2.42, -2.80, -3.25, -3.82,-5.03])
    
    if z == 7:
        mstar = 10**np.array([7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25,10.75, 11.25])
        y = np.array([-1.63,-2.07,-2.49,-2.96,-3.47,-4.11,-4.61,-5.24])
        y_upper = 10**np.array(y + np.array([0.54,0.45,0.38,0.32,0.32,0.41,0.72,0.90]))
        y_lower = 10**np.array(y - np.array([0.54,0.41,0.32,0.30,0.35,0.57,0.82,0.57]))


    if z == 8:
        mstar = 10**np.array([7.25, 7.75, 8.25, 8.75, 9.25, 9.75])
        y = 10**np.array([-1.73,-2.28,-2.88,-3.45,-5.31])
        y_upper = 10**np.array(y + np.array([1.01,0.84,0.75,0.57,0.63,1.01]))
        y_lower = 10**np.array(y - np.array([0.84,0.64,0.47,0.60,0.78,1.64]))
    
    
    ax.fill_between(np.log10(mstar), y_upper,  y_lower,
                         label='Song2016', color=color, alpha=0.5)
    ax.plot(np.log10(mstar), y, color="black")
    print("I did something")
 
    return None

def SMF_Duncan(ax, z, color='magenta'):

    def Schechter(M, Phi, Mstar, alpha):
        return Phi*np.log(10)*10**((M-Mstar)*(alpha+1))*np.exp(-10**((M-Mstar)))

    mstar = np.array([8., 8.5, 9., 9.5, 9.75])
    
    if z not in [4,5,6,7]:
        return None

    if z == 4:
        y = Schechter(mstar, 1.89e-4, 10.51, -1.89)
        y_upper = Schechter(mstar, 5.35e-4, 10.87, -1.74)
        y_lower = Schechter(mstar, 5.7e-5, 10.19, -2.01)
    
    if z == 5:
        y = Schechter(mstar, 2.21e-4, 10.51, -1.64)
        y_upper = Schechter(mstar, 3.01e-4, 10.51, -1.64)
        y_lower = Schechter(mstar, 1.54e-4, 10.51, -1.64)

    if z == 6:
        y = Schechter(mstar, 0.46e-4, 10.51,-1.90)
        y_upper = Schechter(mstar, 8.2e-5, 10.51, -1.90)
        y_lower = Schechter(mstar, 2e-5, 10.51, -1.90)

    if z == 7:
        y = Schechter(mstar, 0.36, 10.51, -1.89)
        y_upper = Schechter(mstar, 3.37e-4, 10.51, -1.89)
        y_lower = Schechter(mstar, 0.01e-4, 10.51, -1.89)
    
    ax.fill_between(np.log10(10**mstar), y_upper, y_lower,
            label='Duncan2014', color=color, alpha=0.3)
    ax.plot(np.log10(10**mstar), y, color='black')
    return None

def SMF_Moustakas(ax, z, color = "blue"):
    if z < 0.5:
        log_mstar = np.array([9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12])
        y_val_raw = np.array([-1.899,-1.923,-1.970,-2.031,-2.055,-2.106,-2.144,-2.179,-2.188,-2.2160,-2.2340,-2.2350,-2.2620,-2.2520,-2.2850,-2.3170,-2.3650,-2.4190,-2.5040,-2.6070,-2.7280,-2.8880,-3.1040,-3.3320,-3.6060,-3.953,-4.363,-4.778,-5.255,-5.87,-6.49])

        ax.scatter(log_mstar,10**y_val_raw,marker="+", color = color, label = "Moustakas 2013")
        return None
    else:
        return None

def SFR_Whitacker(ax, z, color='blue'):

    z_vals = [0.5, 1, 1.5, 2]

    if z not in z_vals:
        return None

    if z == 0.5:
        log_mstar = [8.4,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1]
        log_sfr = [-0.60,-0.38,-0.20,0.05,0.23,0.45,0.65,0.78,0.99,1.06,1.09,1.22,1.21,1.27]
        log_sfr_err = [0.50,0.11,0.05,0.03,0.02,0.02,0.02,0.03,0.07,0.05,0.07,0.07,0.08,0.10]

    if z == 1:
        log_mstar = [8.8,9.1,9.3,9.5,9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1,11.3]
        log_sfr = [-0.03,0.17,0.38,0.64,0.81,1.02,1.18,1.35,1.47,1.58,1.69,1.74,1.81]
        log_sfr_err = [0.13,0.08,0.05,0.02,0.02,0.03,0.04,0.04,0.05,0.05,0.08,0.13,0.11]

    if z == 1.5:
        log_mstar = [9.2,9.4,9.7,9.9,10.1,10.3,10.5,10.7,10.9,11.1,11.3,11.5]
        log_sfr = [0.48,0.70,0.94,1.15,1.38,1.54,1.70,1.83,1.90,2.02,2.19,2.25]
        log_sfr_err = [0.07,0.03,0.02,0.03,0.03,0.04,0.10,0.10,0.11,0.10,0.09,0.18]

    if z == 2:
        log_mstar = [9.3,9.6,9.8,10,10.3,10.5,10.7,10.9,11.1,11.3,11.5]
        log_sfr = [0.82,1.05,1.26,1.46,1.64,1.86,1.95,2.07,2.20,2.32,2.39]
        log_sfr_err = [0.06,0.03,0.03,0.03,0.03,0.05,0.08,0.06,0.10,0.12,0.17]
        

    print(len(log_sfr))
    print(len(log_sfr_err))
    ax.errorbar(log_mstar, log_sfr, yerr=log_sfr_err,color=color, label = "Whitacker 2014")
    ax.legend()

    return None