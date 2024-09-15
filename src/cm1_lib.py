# ----------------------------------------------------------------------------------------------------------------------------
#                  CM1 Library
# Author : Morten Kretschmer
# 
# Read sounding from CM1 input sounding files
# ----------------------------------------------------------------------------------------------------------------------------

from scipy import interpolate
import numpy as np
import csv
from src.meteolib import cr

# ----------------------------------------------------------------------------------------------------------------------------

def read_modelsounding(filename='input_sounding'):
    data_ctr = 0
    Z = []
    Th = []
    r_v = []
    u = []
    v = []
    with open(filename, newline = '') as int_snd:
        snd_dat = csv.reader(int_snd)
        for dataon in snd_dat:
            data_ctr = data_ctr + 1
            if data_ctr == 1: # get surface data from first line
                p_sfc = float(dataon[1])*100.0 
                Th_sfc = float(dataon[2])
                r_sfc = float(dataon[3])/1000.0
                Z.append(0.0)
                Th.append(Th_sfc)
                r_v.append(r_sfc)
                u.append(0.0)
                v.append(0.0)
            else:
                Z.append(float(dataon[0]))
                Th.append(float(dataon[2]))
                r_v.append(float(dataon[3])/1000.0)
                u.append(float(dataon[4]))
                v.append(float(dataon[5]))
    
    return Z, Th, r_v, u, v, p_sfc


# ----------------------------------------------------------------------------------------------------------------------------
# interpolate sounding data onto a regular grid

def interpolating_sounding(Z, Th, r_v, u, v, dz = 100, umove=0.0, vmove=0.0):

    znew = np.arange(0,20*1000, dz, dtype=float)

    f = interpolate.interp1d(Z,Th,fill_value="extrapolate",kind="linear"); Th = f(znew)
    f = interpolate.interp1d(Z,r_v,fill_value="extrapolate",kind="linear"); r_v = f(znew)
    f = interpolate.interp1d(Z,u,fill_value="extrapolate",kind="linear"); u = f(znew) + umove
    f = interpolate.interp1d(Z,v,fill_value="extrapolate",kind="linear"); v = f(znew) + vmove

    return (np.array(znew,dtype=float), np.array(Th,dtype=float), np.array(r_v,dtype=float), 
            np.array(u,dtype=float), np.array(v,dtype=float))


# ----------------------------------------------------------------------------------------------------------------------------
# get pressure using the hydrostatic equation

def calculate_density_potential_temperature(Th, r_v):
    return Th*(1 + r_v/cr['eps'])/(1 + r_v)


def calculate_PII(Z, Th_rho, p_sfc=1000.0):

    Pii = np.zeros(Z.shape[0])
    Pii[0] = (p_sfc/(1000.0*100.0))**(cr['ROCP'])
    intarg = -(cr['G']/cr['cpd'])/Th_rho

    for iz in np.arange(1,Z.shape[0],1):
        Pii[iz] = Pii[iz-1] + 0.5*(intarg[iz]+intarg[iz-1])*(Z[iz]-Z[iz-1])

    return Pii


def calculate_pressure(Pii, p_ref=1000.0):
    return p_ref*100.0*Pii**(cr['cpd']/cr['Rd'])

def calculate_temperature(pres, Th):
    return Th*(pres/(1000*100.0))**(cr['ROCP'])

def calculate_temperature_density(Pii, Th, r_v):
    return Th*Pii*(1 + r_v/cr['eps'])/(1 + r_v)

def calculate_density(pres, T_rho):
    return pres/(cr['Rd']*T_rho)

# ----------------------------------------------------------------------------------------------------------------------------
