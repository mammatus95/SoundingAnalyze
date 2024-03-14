#!/usr/bin/python3

import matplotlib
matplotlib.use('Agg')
# Force matplotlib to not use any Xwindows backend.
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter,ScalarFormatter)
from io import StringIO

import numpy as np
import numpy.ma as ma
import sys
#import datetime
import time

### Constants Used ####
c1 = 0.0498646455 ; c2 = 2.4082965 ; c3 = 7.07475
c4 = 38.9114 ; c5 = 0.0915 ; c6 = 1.2035
eps = 0.62197

MISSING = -9999.0       # Missing Flag
ROCP = 0.28571426       # R over Cp
ZEROCNK = 273.15        # Zero Celsius in Kelvins
G = 9.80665             # Gravity
TOL = 1e-10             # Floating Point Tolerance

##################################################

### functions ###
#Quelle: sharppy

def virtuelle (q, t):
    """
    Returns the virtual temperature in K
    Parameters
    ----------
    w : number, numpy array
        specific humidity (g/kg)
    t : number, numpy array
        temperature (K)
    """
    return t*(1+0.608*q)

def q_to_mixrat(q):
    """
    Parameters
    ----------
    q : number, numpy array
        Specific Humidity (kg/kg)
    """
    return q/(1.0-q)

def mixrat_to_q(mr):
    """
    Parameters
    ----------
    q : number, numpy array
        Mixing Ratio (kg/kg)
    """
    return mr/(1.0+mr)

def temp_at_mixrat(w, p):
    '''
    Returns the temperature (C) of air at the given mixing ratio (g/kg) and
    pressure (hPa)
    Parameters
    ----------
    w : number, numpy array
        Mixing Ratio (g/kg)
    p : number, numpy array
        Pressure (hPa)
    Returns
    -------
    Temperature (C) of air at given mixing ratio and pressure
    '''
    x = np.log10(w * p / (622. + w))
    x = (np.power(10.,((c1 * x) + c2)) - c3 + (c4 * np.power((np.power(10,(c5 * x)) - c6),2))) - ZEROCNK
    return x

def satdampfdruck (t):
    #wahrscheindlich der Saettigungsdampfdruck fur wasser also der maximale partialdruck den der wasserdampf erreichen kann
    pol = t * (1.1112018e-17 + (t * -3.0994571e-20))
    pol = t * (2.1874425e-13 + (t * (-1.789232e-15 + pol)))
    pol = t * (4.3884180e-09 + (t * (-2.988388e-11 + pol)))
    pol = t * (7.8736169e-05 + (t * (-6.111796e-07 + pol)))
    pol = 0.99999683 + (t * (-9.082695e-03 + pol))
    return 6.1078 / pol**8

def saettigungsdampfdruck (t):
     a = 7.5
     b = 237.3
     return 6.1078 * 10**((a*t)/(b+t))

def wobf(t):

    t = t - 20
    if type(t) == type(np.array([])) or type(t) == type(np.ma.array([])):
        npol = 1. + t * (-8.841660499999999e-3 + t * ( 1.4714143e-4 + t * (-9.671989000000001e-7 + t * (-3.2607217e-8 + t * (-3.8598073e-10)))))
        npol = 15.13 / (np.power(npol,4))
        ppol = t * (4.9618922e-07 + t * (-6.1059365e-09 + t * (3.9401551e-11 + t * (-1.2588129e-13 + t * (1.6688280e-16)))))
        ppol = 1 + t * (3.6182989e-03 + t * (-1.3603273e-05 + ppol))
        ppol = (29.93 / np.power(ppol,4)) + (0.96 * t) - 14.8
        correction = np.zeros(t.shape, dtype=np.float64)
        correction[t <= 0] = npol[t <= 0]
        correction[t > 0] = ppol[t > 0]
        return correction
    else:
        if t is np.ma.masked:
            return t
        if t <= 0:
            npol = 1. + t * (-8.841660499999999e-3 + t * ( 1.4714143e-4 + t * (-9.671989000000001e-7 + t * (-3.2607217e-8 + t * (-3.8598073e-10)))))
            npol = 15.13 / (np.power(npol,4))
            return npol
        else:
            ppol = t * (4.9618922e-07 + t * (-6.1059365e-09 + t * (3.9401551e-11 + t * (-1.2588129e-13 + t * (1.6688280e-16)))))
            ppol = 1 + t * (3.6182989e-03 + t * (-1.3603273e-05 + ppol))
            ppol = (29.93 / np.power(ppol,4)) + (0.96 * t) - 14.8
            return ppol

def satlift(p, thetam):

    #if type(p) == type(np.array([p])) or type(thetam) == type(np.array([thetam])):
    if np.fabs(p - 1000.) - 0.001 <= 0: return thetam
    eor = 999
    while np.fabs(eor) - 0.1 > 0:
        if eor == 999:                  # First Pass
            pwrp = np.power((p / 1000.),ROCP)
            t1 = (thetam + ZEROCNK) * pwrp - ZEROCNK
            e1 = wobf(t1) - wobf(thetam)
            rate = 1
        else:                           # Successive Passes
            rate = (t2 - t1) / (e2 - e1)
            t1 = t2
            e1 = e2
        t2 = t1 - (e1 * rate)
        e2 = (t2 + ZEROCNK) / pwrp - ZEROCNK
        e2 += wobf(t2) - wobf(e2) - thetam
        eor = e2 * rate
    return t2 - eor

def lcltemp(t, td):
    '''
    Returns the temperature (C) of a parcel when raised to its LCL.
    Parameters
    ----------
    t : number, numpy array
        Temperature of the parcel (C)
    td : number, numpy array
        Dewpoint temperature of the parcel (C)
    Returns
    -------
    Temperature (C) of the parcel at it's LCL.
    '''
    s = t - td
    dlt = s * (1.2185 + 0.001278 * t + s * (-0.00219 + 1.173e-5 * s -
        0.0000052 * t))
    return t - dlt


def thalvl(theta, t):
    '''
    Returns the level (hPa) of a parcel.
    Parameters
    ----------
    theta : number, numpy array
        Potential temperature of the parcel (C)
    t : number, numpy array
        Temperature of the parcel (C)
    Returns
    -------
    Pressure Level (hPa [float]) of the parcel
    '''
    t = t + ZEROCNK
    theta = theta + ZEROCNK
    return 1000. / (np.power((theta / t),(1./ROCP)))


def theta(p, t, p2=1000.):
    '''
    Returns the potential temperature (C) of a parcel.
    Parameters
    ----------
    p : number, numpy array
        The pressure of the parcel (hPa)
    t : number, numpy array
        Temperature of the parcel (C)
    p2 : number, numpy array (default 1000.)
        Reference pressure level (hPa)
    Returns
    -------
    Potential temperature (C)
    '''
    return ((t + ZEROCNK) * np.power((p2 / p),ROCP)) - ZEROCNK

def thetas(theta, presvals):
    return ((theta + ZEROCNK) / (np.power((1000. / presvals),ROCP))) - ZEROCNK

def drylift(p, t, td):
    '''
    Lifts a parcel to the LCL and returns its new level and temperature.
    Parameters
    ----------
    p : number, numpy array
        Pressure of initial parcel in hPa
    t : number, numpy array
        Temperature of inital parcel in C
    td : number, numpy array
        Dew Point of initial parcel in C
    Returns
    -------
    p2 : number, numpy array
        LCL pressure in hPa
    t2 : number, numpy array
        LCL Temperature in C
    '''
    t2 = lcltemp(t, td)
    p2 = thalvl(theta(p, t, 1000.), t2)
    return p2, t2

def wetlift(p, t, p2):
    thta = theta(p, t, 1000.)
    thetam = thta - wobf(thta) + wobf(t)
    return satlift(p2, thetam)

def wetbulb(p, t, td):
    '''
    Calculates the wetbulb temperature (C) for the given parcel
    Parameters
    ----------
    p : number
        Pressure of parcel (hPa)
    t : number
        Temperature of parcel (C)
    td : number
        Dew Point of parcel (C)
    Returns
    -------
    Wetbulb temperature (C)
    '''
    p2, t2 = drylift(p, t, td)
    return wetlift(p2, t2, p)

b, ax = plt.subplots(1,figsize=(15,11))

ax.grid(True)

pmax = 1000
pmin = 100
dp = -10
presvals = np.arange(pmax, pmin+dp, dp)

# plot the moist-adiabats
for t in np.array([-40, -30, -20, -10,  0,  5, 10, 15, 20, 25, 30, 35, 40]):
    tw = []
    for p in presvals:
        tw.append(wetlift(1000., t, p))
    ax.semilogy(tw, presvals, 'b--',lw=0.5, alpha=0.7)

# plot the dry adiabats
for t in np.array([-80, -60, -40, -30, -20, -10,   0,  10,  20,  30,  40, 50,  70,  90, 110, 130, 150, 170, 190, 210, 230]):
    ax.semilogy(thetas(t, presvals), presvals, 'r--',lw=0.5, alpha=0.7)

# linie gleichensaettiungsmischungsverhaeltnis
pmin = 400
presvals = np.arange(pmax, pmin+dp, dp)

mrr = np.array([0.1,0.5,1.0,5.0,10.0,15.0,20.0,30.0])
for i in range(0,np.size(mrr)):
    ax.semilogy(temp_at_mixrat(mrr[i], presvals),presvals,'m--',lw=0.6, alpha=0.8)

plt.yscale('log')

ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_yticks(np.linspace(100,1000,10))
ax.set_ylim(1050,400)
ax.set_xticks(np.arange(-90,50,10))
ax.set_xlim(-90,50)

ax.set_title("Stuve Diagramm",fontsize=18)
plt.xlabel("Temperatur [C]")
plt.ylabel("Druck [hPa]")
plt.savefig("stuve2.png")
plt.close()
