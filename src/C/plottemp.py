#!/usr/bin/python3

import matplotlib
matplotlib.use('Agg')
# Force matplotlib to not use any Xwindows backend.
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter,ScalarFormatter)
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from io import StringIO

import numpy as np
import numpy.ma as ma
import sys
#import datetime
import time
from copy import deepcopy

from skewT import *

"""
def datum ():
    today = datetime.date.today()
    today = today.timetuple()
    hour = int(time.strftime("%H"))

    x = datetime.datetime(today.tm_year, today.tm_mon, today.tm_mday, 0)

    uhrzeit = " Uhrzeit: " + x.strftime("%H:%M") + " UTC "
    return x.strftime("%d.%m.%Y ") + uhrzeit
"""
def datum ():
    hour = int(time.strftime("%H"))
    if hour >= 18:
        return time.strftime("%d.%m.%Y Time: 18 UTC")
    elif hour >= 12:
        return time.strftime("%d.%m.%Y Time: 12 UTC")
    elif hour >= 6:
        return time.strftime("%d.%m.%Y Time: 06 UTC")
    else:
        return time.strftime("%d.%m.%Y Time: 00 UTC")



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

def wetpot(p,t,p2=1000.):
    return wetlift(p, t, p2)


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
#Quelle: SHARPpy
def thetae(p,t,q,p2=1000.):
    """
    Calculates the equlivalent potential temperature (K) for the given parcel
    Parameters
    ----------
    p : number
        Pressure of parcel (hPa)
    t : number
        Temperature of parcel (C)
    q : number
        specific humidity (kg/kg)
    Returns
    -------
    qulivalent potential temperature (K)
    """
    return ((t + ZEROCNK)+((2.501*1000000)/1004)*q)* np.power((p2/p),ROCP)

def mean_wind(u,v,ps, stu=0, stv=0):
    return np.average(u, weights=ps)-stu, np.average(v, weights=ps)-stv

def non_parcel_bunkers_motion_experimental(u,v,ps,i_500m,i_5km,i_6km):
    d=7.5
    ## sfc-500m Mean Wind
    mnu500m, mnv500m = mean_wind(u[:i_500m],v[:i_500m],ps[:i_500m])
    
    ## 5.5km-6.0km Mean Wind
    mnu5500m_6000m, mnv5500m_6000m = mean_wind(u[i_5km:i_6km],v[i_5km:i_6km],ps[i_5km:i_6km])
    
    # shear vector of the two mean winds
    shru = mnu5500m_6000m - mnu500m
    shrv = mnv5500m_6000m - mnv500m
    
    # SFC-6km Mean Wind
    mnu6, mnv6 =  mean_wind(u[:i_6km],v[:i_6km],ps[:i_6km])
    
    # Bunkers Right Motion
    tmp = d / np.sqrt(shru*shru + shrv*shrv)
    rstu = mnu6 + (tmp * shrv)
    rstv = mnv6 - (tmp * shru)
    lstu = mnu6 - (tmp * shrv)
    lstv = mnv6 + (tmp * shru)
    
    return rstu, rstv, lstu, lstv, mnu6, mnv6

#################################################################################

###load sounding data

datafile = open('radiodata.txt', 'r').read()

data = np.array([l.strip() for l in datafile.split('\n')])
#height, p_temp, temp, td_temp, mix, virtemp, pot, pseudo_pot, laps, cape, cin, speed, u, v, srur, srvr, ushr, vshr, srul, srvl = np.genfromtxt(data, delimiter=',', unpack=True ) #issue #3184
height, p_temp, t_temp, td_temp, mix, vir_temp, pot, pseudo_pot, laps, cape, cin, speed, u, v, srur, srvr, ushr, vshr, srul, srvl = np.loadtxt(data, delimiter=',', unpack=True )
laps[0] = 0
cin = cin *(-1)
anzahl = np.size(height)
##################################################################################
print(height)
#system arguements
#sys.argv[1] stationsname
#sys.argv[2] ret
#sys.argv[3] mu cape
#sys.argv[4] mu lclp
#sys.argv[5] mu lclt
#sys.argv[6] rstu
#sys.argv[7] rstv
#sys.argv[8] lstu
#sys.argv[9] lstv
#sys.argv[10] efi bot
#sys.argv[11] efi top
#sys.argv[12] wxshear
#sys.argv[13] brn ml
#sys.argv[14] lapse 850
#sys.argv[15] mix mean
print (sys.argv)

if sys.argv[1] == '10035':
    station = 'Schleswig'
elif sys.argv[1] == '10113':
    station = 'Norderney'
elif sys.argv[1] == '10184':
    station = 'Greifswald'
elif sys.argv[1] == '10238':
    station = 'Bergen'
#elif sys.argv[1] == '10304':
#    station = ''
elif sys.argv[1] == '10393':
    station = 'Lindenberg'
elif sys.argv[1] == '10410':
    station = 'Essen'
elif sys.argv[1] == '10548':
    station = ' Meiningen'
elif sys.argv[1] == '10618':
    station = 'Idar-Oberstein'
elif sys.argv[1] == '10771':
    station = 'Kuemmersbruck'
elif sys.argv[1] == '10739':
    station = 'Stuttgart'
elif sys.argv[1] == '10868':
    station = 'Muenchen-Oberschlssheim'
#elif sys.argv[1] == '06610':
#    station = ''
elif sys.argv[1] == '11520':
    station = 'Praha-Libus'
elif sys.argv[1] == '11035':
    station = 'Wien'
elif sys.argv[1] == '16045':
    station = 'Rivolto'
elif sys.argv[1] == '11747':
    station = 'Prostejov'
else:
    station = sys.argv[1]

title = "Sounding of " + station 
wenn = "Date: " + datum()

name1 = sys.argv[1] + "thermo.png"
name2 = sys.argv[1] + "thermo_3000.png"
name3 = sys.argv[1] + "theta.png"
name4 = sys.argv[1] + "wind.png"
name5 = sys.argv[1] + "skewT.png"
name6 = sys.argv[1] + "nixon.png"
name7 = sys.argv[1] + "ri.png"

##################################################################################

#print (pot[int(sys.argv[4])], temp[int(sys.argv[4])])
#temp.mu.cape,drk[ret].pottemp,drk[ret].temp, (int) drk[ret].lclp, drk[ret].lclt, ret, sharppy.rstu,sharppy.rstv,sharppy.rstu,sharppy.rstv

#wetbulb
wetzero = np.zeros(np.size(t_temp))
for i in range(0,np.size(t_temp)):
    #print(p_temp[i],t_temp[i],td_temp[i])
    wetzero[i] = wetbulb(p_temp[i],t_temp[i],td_temp[i])

#virtuell temperature vir_temp
#q_temp = mixrat_to_q(mix/1000.0)
#tv_temp = virtuelle (q_temp, t_temp)

#lclp, theoretical parcel
"""
pmax = int(p_temp[0])
pmin = 100
dp = -2
press_parcel2 = np.arange(pmax, pmin+dp, dp)
lclp, lclt = drylift(p_temp[0], 22.5, 14.8)# sb cape
mu2_parcel = np.zeros(press_parcel2.size, dtype=np.float64)
for i in range(0,press_parcel2.size):
    if press_parcel2[i] >= lclp:
        #mu_parcel[i] = thetas(theta(p_temp[int(sys.argv[2])], t_temp[int(sys.argv[2])], p2=1000.), press_parcel[i])
        mu2_parcel[i] = thetas(theta(lclp, lclt, p2=1000.), press_parcel2[i])
    else:
        mu2_parcel[i] = wetlift(lclp,lclt,press_parcel2[i])

"""
pmax = int(p_temp[int(sys.argv[2])])
pmin = 100
dp = -2
press_parcel = np.arange(pmax, pmin+dp, dp)

lclp = float(sys.argv[4])
lclt = float(sys.argv[5])
#mu = np.argmax(pseudo_pot)
#print(mu,pseudo_pot)
#lclp, lclt = drylift(p_temp[mu], t_temp[mu], td_temp[mu])# sb cape

mu_parcel = np.zeros(press_parcel.size, dtype=np.float64)
for i in range(0,press_parcel.size):
    if press_parcel[i] >= lclp:
        #mu_parcel[i] = thetas(theta(p_temp[int(sys.argv[2])], t_temp[int(sys.argv[2])], p2=1000.), press_parcel[i])
        mu_parcel[i] = thetas(theta(lclp, lclt, p2=1000.), press_parcel[i])
    else:
        mu_parcel[i] = wetlift(lclp,lclt,press_parcel[i])

        
    
#############################################
 
"""      
string = 'LCL: '
string += sys.argv[4]
string += ' hPa'

string += ', MUCAPE: '
string += sys.argv[3]
string += 'J/kg                    '
"""

b, ax = plt.subplots(1,figsize=(15,11))
plt.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.96)
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
    ax.semilogy(tw, presvals, 'b--',lw=0.7, alpha=0.7)

# plot the dry adiabats
for t in np.array([-80, -60, -40, -30, -20, -10,   0,  10,  20,  30,  40, 50,  70,  90, 110, 130, 150, 170, 190, 210, 230]):
    ax.semilogy(thetas(t, presvals), presvals, 'r--',lw=0.7, alpha=0.7)

# linie gleichensaettiungsmischungsverhaeltnis
pmin = 700
presvals = np.arange(pmax, pmin+dp, dp)

mrr = np.append(mix[0], np.array([0.1,0.5,1.0,5.0,10.0,15.0,20.0,30.0]))
for i in range(0,np.size(mrr)):
    ax.semilogy(temp_at_mixrat(mrr[i], presvals),presvals,'m--',lw=0.6, alpha=0.8)

#plot souding datas
ax.plot(wetzero,p_temp,'-c',lw=1.1)
ax.plot(vir_temp,p_temp,'-b',lw=1.1)
ax.plot(t_temp,p_temp,'-b',td_temp,p_temp,'-g',lw=1.5)
ax.plot(mu_parcel,press_parcel,'-r',lw=1.8)
#ax.plot(mu2_parcel,press_parcel2,'-m',lw=1.9)

if float(sys.argv[3]) > 10.0:
    mu = np.zeros(p_temp.size, dtype=np.float64)
    for i in range(int(sys.argv[2]),p_temp.size):
        #if p_temp[i] >= lclp:
        #    mu[i] = thetas(theta(p_temp[int(sys.argv[2])], t_temp[int(sys.argv[2])], p2=1000.), p_temp[i])
        #else:
        mu[i] = wetlift(lclp,lclt,p_temp[i])

    ax.fill_betweenx(p_temp,t_temp, mu,  where=mu >= t_temp, facecolor='red', alpha=0.4)


plt.yscale('log')

ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_yticks(np.linspace(100,1000,10))
ax.set_ylim(1050,100)
ax.set_xticks(np.arange(-90,50,10))
ax.set_xlim(-90,50)

#plt.annotate(string, xy=(float(sys.argv[5])+20, float(sys.argv[4])), xytext=(float(sys.argv[5]) - 85, float(sys.argv[4])+10), arrowprops=dict(arrowstyle="-",facecolor='black'))
ax.axhline(float(sys.argv[10]))#efi bot
ax.axhline(float(sys.argv[11]))#efi top

#wind bars
anzahl_bar=45
anzahl=(np.where(p_temp==100))[0][0]
iter_bar = int(anzahl/anzahl_bar)
x = np.zeros(np.size(p_temp[:anzahl:iter_bar])) + 45
ax.barbs(x, p_temp[:anzahl:iter_bar], u[:anzahl:iter_bar], v[:anzahl:iter_bar], np.sqrt(u[:anzahl:iter_bar] * u[:anzahl:iter_bar] + v[:anzahl:iter_bar] * v[:anzahl:iter_bar]), length=8, pivot='middle',barb_increments=dict(half=2.5, full=5, flag=25))

plt.text(-90,99,wenn,fontsize=11)
ax.set_title(title,fontsize=18)
plt.xlabel("Temperatur [C]")
plt.ylabel("Druck [hPa]")
plt.savefig(name1)
plt.close()

#############################################

"""
f, axarr = plt.subplots(2,2,figsize=(8,8))
axarr[0,0].plot(cape,height,'-r',cin,height,'-g',cape,height,'-m',lw=1.5)
axarr[0,0].set_title('CAPE/CIN in J/kg')
max_cape = float(sys.argv[3])+100
if max_cape < 100:
    max_cape = 100
axarr[0,0].axis([0,max_cape,0,3000])
axarr[0,0].grid()

axarr[0,1].plot((laps[2:]+laps[:-2])/2.0,height[1:-1],'-g',lw=1.2)
axarr[0,1].set_title('Lapse Rates in K/KM')
axarr[0,1].axvline(9.8, color='k')
axarr[0,1].axis([-30,30,0,3000])
axarr[0,1].grid()

axarr[1,0].grid()

pmax = 1050
pmin = 650
dp = -10
presvals = np.arange(pmax, pmin+dp, dp)

# plot the moist-adiabats
for t in np.array([-40, -30, -20, -10,  0,  5, 10, 15, 20, 25, 30, 35, 40]):
    tw = []
    for p in presvals:
        tw.append(wetlift(1000., t, p))
    axarr[1,0].semilogy(tw, presvals, 'b--',lw=0.5, alpha=0.7)

# plot the dry adiabats
for t in np.array([-80, -60, -40, -30, -20, -10,   0,  10,  20,  30,  40, 50,  70,  90, 110, 130, 150, 170, 190, 210, 230]):
    axarr[1,0].semilogy(thetas(t, presvals), presvals, 'r--',lw=0.5, alpha=0.7)

# linie gleichensaettiungsmischungsverhaeltnis
mrr = np.append(mix[0], np.array([0.1,0.5,1.0,5.0,10.0,15.0,20.0,30.0]))
for i in range(0,np.size(mrr)):
    axarr[1,0].semilogy(temp_at_mixrat(mrr[i], presvals),presvals,'m--',lw=0.6, alpha=0.8)

axarr[1,0].semilogy(t_temp,p_temp,'-b',td_temp,p_temp,'--b',lw=1.3)
axarr[1,0].semilogy(wetzero,p_temp,'-c',lw=0.9)
axarr[1,0].semilogy(vir_temp,p_temp,'-g',lw=0.9)
axarr[1,0].semilogy(mu_parcel,press_parcel,'-r',lw=1.5)
axarr[1,0].set_title('Temperatur/Feuchteprofil')

axarr[1,0].axis([t_temp[0]-50,t_temp[0] + 25,1030,650])
x = np.zeros(np.size(p_temp[::5])) + (t_temp[0] + 20)
axarr[1,0].barbs(x, p_temp[::5], u[::5], v[::5], np.sqrt(u[::5] * u[::5] + v[::5] * v[::5]), length=6, pivot='middle')
plt.yscale('log')

axarr[1,1].semilogy(mix,p_temp,'-r',lw=1.5)
axarr[1,1].set_title('Mixing ratio g/kg')
axarr[1,1].axis([-10,35,1030,650])
axarr[1,1].set_xticks(np.arange(0,35,5))
axarr[1,1].grid()
plt.yscale('log')
plt.savefig(name2)

#plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
"""

#############################################

pote = thetae(p_temp,t_temp,mixrat_to_q(mix/1000.0))-ZEROCNK

b, ax = plt.subplots(1)

ax.plot(pot, p_temp,'-b', pseudo_pot, p_temp,'-r', wetpot(p_temp,t_temp), p_temp, '-g',pote, p_temp,'-m',lw=1.1)
plt.legend(('potential temperature ','pseudo pot. temperature','wetbulb pot. temperature','equivalent pot. temperature' ),loc=4)
ax.set_title(title)
plt.yscale('log')
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_yticks(np.linspace(100,1000,10))
ax.set_ylim(1050,100)
ax.set_xticks(np.append(np.arange(-90,90,10,dtype=np.int64),np.arange(100,300,20,dtype=np.int64)))
ax.set_xlim(pot[0] - 15,pot[0] + 180)

plt.xlabel("Temperatur [C]")
plt.ylabel("Druck [hPa]")
plt.grid(True)

#ax2 = ax.twiny()
#ax2.set_xlabel("Temperatur [K]")
#ax2.set_xlim(pot[0] + ZEROCNK, 200 + ZEROCNK)
#ax2.set_xticks(np.arange(pot[0] - 10 + ZEROCNK,200 + ZEROCNK,50,dtype=np.int64))
#ax2.set_xticks(np.arange(260,500,40,dtype=np.int64))
#ax2.set_xticklabels(['7','8','99'])


plt.savefig(name3)
plt.close()

#Kelvin-Helmholz Plot

def df_dz(fmat,zh):
    """
    Parameter:
    zh  : major-layer height of the Modell
    fmat: field of which we want the derivation
    """
    df_dz = np.zeros((fmat.shape))
    df_dz[1:-1] = (fmat[:-2]-fmat[2:])/(zh[:-2]-zh[2:])
    #Rand:
    df_dz[0]    = (fmat[1]-fmat[0])/(zh[1]-zh[0])
    df_dz[-1]   = (fmat[-1]-fmat[-2])/(zh[-1]-zh[-2])
    return df_dz


pot_z=df_dz(pot,height)
pote_z=df_dz(pote,height)
Ri = ((G/pot)*pot_z)/(df_dz(u,height)**2 + df_dz(v,height)**2)

fig, axs = plt.subplots(1, 4, sharey=True, tight_layout=True)
plt.yscale('log')
axs[0].plot(Ri,p_temp)
axs[0].axvline(x=0,c='k')
axs[1].plot(pot_z,p_temp)
axs[1].axvline(x=0,c='k')
axs[2].plot(df_dz(u,height),p_temp)
axs[2].axvline(x=0,c='k')
axs[3].plot(G*((pote[0]-pot)/(pot+ZEROCNK)),p_temp,'r-')
axs[3].plot(G*((pseudo_pot[0]-pot)/(pot+ZEROCNK)),p_temp,'m-')
axs[3].axvline(x=0,c='k')
axs[3].set_xlim(-0.5,1)
#axs[0].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axs[0].set_yticks(np.linspace(100,1000,10))
axs[0].set_ylim(1050,300)
plt.savefig(name7)
plt.close()

#############################################
#                                           #
#              Hodograph                    #
#                                           #
#############################################

b, ax4 = plt.subplots(212)
ax4 = plt.axes()
ax4.get_xaxis().set_visible(False)
ax4.get_yaxis().set_visible(False)
for i in range(10,100,10):
    # Draw the range rings around the hodograph.
    circle = plt.Circle((0,0),i,color='k',alpha=.3, fill=False)
    ax4.add_artist(circle)

ax4.set_xlim(-50,50)
ax4.set_ylim(-50,50)
ax4.set_title(title)
ax4.axhline(y=0, color='k')
ax4.axvline(x=0, color='k')

ax4.text(-2.45,40.2,"40 m/s")
ax4.text(-2.45,20.2,"20 m/s")

ax4.text(-2.45,-39.8,"40 m/s")
ax4.text(-2.45,-19.8,"20 m/s")

m=[]
n=[]
i=0
while height[i] < 500 :
    m.append(u[i])
    n.append(v[i])
    i += 1
i_500m = deepcopy(i)
while height[i] < 1000 :
	m.append(u[i])
	n.append(v[i])
	i += 1
u_1 = u[i]
v_1 = v[i]
m.append(u[i])
n.append(v[i])
#print("%.1f %.1f" % (m,n))
ax4.plot(m,n, 'g-', lw=2, label='0-1 km')
#ax4.plot(m,n, 'gx', lw=2)
m=[]
n=[]
while height[i] < 3000 :
	m.append(u[i])
	n.append(v[i])
	i += 1

m.append(u[i])
n.append(v[i])
#print("%.1f %.1f" % (m,n))
ax4.plot(m,n, 'r-', lw=2, label='1-3 km')
#ax4.plot(m[::3],n[::3], 'rx', lw=2)
m=[]
n=[]
while height[i] < 5500:
    m.append(u[i])
    n.append(v[i])
    i += 1
i_5km = deepcopy(i)
while height[i] < 6000 :
	m.append(u[i])
	n.append(v[i])
	i += 1
i_6km = deepcopy(i)
m.append(u[i])
n.append(v[i])
u_6 = u[i]
v_6 = v[i]
#print("%.1f %.1f" % (m,n))
ax4.plot(m,n, 'b-', lw=2, label='3-6 km')
#ax4.plot(m[::3],n[::3], 'bx', lw=2)

m=[]
n=[]

while height[i] < 12000 :
        m.append(u[i])
        n.append(v[i])
        i += 1
m.append(u[i])
n.append(v[i])
#print("%.1f %.1f" % (m,n))
ax4.plot(m,n, 'm-', lw=2, label='6-12 km')
#ax4.plot(m[::3],n[::3], 'bx', lw=2)

ax4.plot(float(sys.argv[6]), float(sys.argv[7]), 'ro') # Plot Bunker's Storm motion right mover as a red dot
ax4.plot(float(sys.argv[8]), float(sys.argv[9]), 'bo') # Plot Bunker's Storm motion left mover as a blue dot

rstu, rstv, lstu, lstv, mwu6, mwv6 = non_parcel_bunkers_motion_experimental(u,v,p_temp,i_500m,i_5km,i_6km)
ax4.plot(mwu6, mwv6, 'o', color = 'grey',alpha=0.7,zorder=2) # 0-6km Mean weight wind

ax4.legend(loc=3)
plt.savefig(name4)
plt.close()


#############################################################
#nixon strom relative

def nixon_hodo(u,v,rstu,rstv):
    #float(sys.argv[6]), float(sys.argv[7])
    """
    u,v : horizontal wind components
    rstu,rstv : storm motion vector for right mover
    """
    return u-rstu,v-rstv

u,v = nixon_hodo(u,v,float(sys.argv[6]), float(sys.argv[7]))

b, ax4 = plt.subplots(212)
ax4 = plt.axes()
ax4.get_xaxis().set_visible(False)
ax4.get_yaxis().set_visible(False)
for i in range(10,100,10):
    # Draw the range rings around the hodograph.
    circle = plt.Circle((0,0),i,color='k',alpha=.3, fill=False)
    ax4.add_artist(circle)

ax4.set_xlim(-50,50)
ax4.set_ylim(-50,50)
ax4.set_title(title)
ax4.axhline(y=0, color='k')
ax4.axvline(x=0, color='k')

ax4.text(-2.45,40.2,"40 m/s")
ax4.text(-2.45,20.2,"20 m/s")

ax4.text(-2.45,-39.8,"40 m/s")
ax4.text(-2.45,-19.8,"20 m/s")

m=[]
n=[]
i=0
while height[i] < 500 :
    m.append(u[i])
    n.append(v[i])
    i += 1
i_500m = deepcopy(i)
while height[i] < 1000 :
	m.append(u[i])
	n.append(v[i])
	i += 1
u_1 = u[i]
v_1 = v[i]
m.append(u[i])
n.append(v[i])
#print("%.1f %.1f" % (m,n))
ax4.plot(m,n, 'g-', lw=2, label='0-1 km')
#ax4.plot(m,n, 'gx', lw=2)
m=[]
n=[]
while height[i] < 3000 :
	m.append(u[i])
	n.append(v[i])
	i += 1

m.append(u[i])
n.append(v[i])
#print("%.1f %.1f" % (m,n))
ax4.plot(m,n, 'r-', lw=2, label='1-3 km')
#ax4.plot(m[::3],n[::3], 'rx', lw=2)
m=[]
n=[]
while height[i] < 5500:
    m.append(u[i])
    n.append(v[i])
    i += 1
i_5km = deepcopy(i)
while height[i] < 6000 :
	m.append(u[i])
	n.append(v[i])
	i += 1
i_6km = deepcopy(i)
m.append(u[i])
n.append(v[i])
u_6 = u[i]
v_6 = v[i]
#print("%.1f %.1f" % (m,n))
ax4.plot(m,n, 'b-', lw=2, label='3-6 km')
#ax4.plot(m[::3],n[::3], 'bx', lw=2)

m=[]
n=[]

while height[i] < 12000 :
        m.append(u[i])
        n.append(v[i])
        i += 1
m.append(u[i])
n.append(v[i])
#print("%.1f %.1f" % (m,n))
ax4.plot(m,n, 'm-', lw=2, label='6-12 km')
#ax4.plot(m[::3],n[::3], 'bx', lw=2)

ax4.arrow(0,0,float(sys.argv[6]), float(sys.argv[7]), head_width=0.05, head_length=0.1, color='black')

#ax4.plot(float(sys.argv[6]), float(sys.argv[7]), 'ro') # Plot Bunker's Storm motion right mover as a red dot
#ax4.plot(float(sys.argv[8]), float(sys.argv[9]), 'bo') # Plot Bunker's Storm motion left mover as a blue dot

#rstu, rstv, lstu, lstv, mwu6, mwv6 = non_parcel_bunkers_motion_experimental(u,v,p_temp,i_500m,i_5km,i_6km)
#ax4.plot(mwu6, mwv6, 'o', color = 'grey',alpha=0.7,zorder=2) # 0-6km Mean weight wind

ax4.legend(loc=3)
plt.savefig(name6)
plt.close()




#############################################################

"""
### SkewT-logP ####
register_projection(SkewXAxes)
#factor=1
#fig = plt.figure(figsize=(6.5875*factor, 6.2125*factor))
fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111, projection='skewx')
#ax = fig.add_subplot(gs[0,0], projection='skewx')
#plt.subplots_adjust(left=0.1,right=0.6)

ax.grid(True)

pmax = 1000
pmin = 10
dp = -10
presvals = np.arange(int(pmax), int(pmin)+dp, dp)

# plot the moist-adiabats
for t in np.arange(-10,45,5):
    tw = []
    for p in presvals:
        tw.append(wetlift(1000., t, p))
    ax.semilogy(tw, presvals, 'k-', alpha=.2)

def thetas(theta, presvals):
    return ((theta + ZEROCNK) / (np.power((1000. / presvals),ROCP))) - ZEROCNK

# plot the dry adiabats
for t in np.arange(-50,110,10):
    ax.semilogy(thetas(t, presvals), presvals, 'r-', alpha=.2)

plt.title(title,fontsize=11)
# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dicatated by the typical meteorological plot
ax.semilogy(wetzero,p_temp,'c-',lw=1.1) # Plot the wetbulb profile
ax.semilogy(vir_temp,p_temp,'-b',lw=1.1) # plot the dewpoint profile
ax.semilogy(t_temp,p_temp,'-b',td_temp,p_temp,'-g',lw=1.5) # Plot the temperature profile
ax.semilogy(mu_parcel,press_parcel,'-r',lw=1.8) # plot the parcel trace
#ax.plot(mu2_parcel,press_parcel2,'-m',lw=1.9)

if float(sys.argv[3]) > 10.0:
    ax.fill_betweenx(p_temp,t_temp, mu,  where=mu >= t_temp, facecolor='red', alpha=0.4)#interpolate=True,

# An example of a slanted line at constant X
l = ax.axvline(0, color='b', linestyle='--')
l = ax.axvline(-20, color='b', linestyle='--')

# Plot the effective inflow layer using blue horizontal lines
ax.axhline(float(sys.argv[10]), color='b')
ax.axhline(float(sys.argv[11]), color='b')

# Disables the log-formatting that comes with semilogy
ax.yaxis.set_major_formatter(plt.ScalarFormatter())
ax.set_yticks(np.linspace(100,1000,10))
ax.set_ylim(1050,100)
ax.xaxis.set_major_locator(plt.MultipleLocator(10))
ax.set_xlim(-40,30)

# Draw the hodograph on the Skew-T.
# TAS 2015-4-16: hodograph doesn't plot for some reason ...
ax2 = plt.axes([.625,.40,.25,.50])

ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
for i in range(10,90,10):
    # Draw the range rings around the hodograph.
    circle = plt.Circle((0,0),i,color='k',alpha=.3, fill=False)
    ax2.add_artist(circle)

i=0
while height[i] < 10000 :
        i += 1
points = np.array([u[:i+1], v[:i+1]]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
# Create a continuous norm to map from data points to colors
#cmap = ListedColormap(['r', 'g', 'b'])
norm = plt.Normalize(0, 10000)
lc = LineCollection(segments, cmap='hsv_r', norm=norm)
# Set the values used for colormapping
lc.set_array(height)
lc.set_linewidth(2)
line = ax2.add_collection(lc)
cbar=fig.colorbar(line, ax=ax2, orientation='horizontal', ticks=np.arange(0,10000,1000))
clevs=np.arange(0,10,1)
cbar.ax.set_xticklabels(labels=clevs[::1])

ax2.plot(float(sys.argv[6]), float(sys.argv[7]), 'ro', lw=0.7,zorder=3) # Plot Bunker's Storm motion right mover as a red dot
ax2.plot(float(sys.argv[8]), float(sys.argv[9]), 'bo', lw=0.7,zorder=3) # Plot Bunker's Storm motion left mover as a blue dot

ax2.set_xlim(-40,40)
ax2.set_ylim(-40,40)
ax2.axhline(y=0, color='k', lw=0.9,zorder=2)
ax2.axvline(x=0, color='k', lw=0.9,zorder=2)
fig.text(0.750, 0.38, "km")
fig.text(0.625, 0.33, "MU CAPE : " + sys.argv[3] + " J/kg")
fig.text(0.625, 0.30, "Mean MixRatio : "+ sys.argv[15] + " g/kg")
fig.text(0.625, 0.27, "BRN ML : " + sys.argv[13])
fig.text(0.625, 0.24, "Bulk Shear 0-1km : %.1f m/s" % ( np.sqrt((u_1-u[0])**2 + (v_1-v[0])**2) ))
fig.text(0.625, 0.21, "Bulk Shear 0-6km : %.1f m/s" % ( np.sqrt((u_6-u[0])**2 + (v_6-v[0])**2) ))
fig.text(0.625, 0.18, "WMAXSHEAR : " + sys.argv[12] + r" $m^2/s^2$")
fig.text(0.625, 0.15, "Lapse Rate (850-500hPa) : "+ sys.argv[14] + " K/km")
#sys.argv[12] wxshear
#sys.argv[13] brn ml
#sys.argv[14] lapse 850
#sys.argv[15] mix mean

plt.savefig(name5)
plt.close()
"""

# Create a SkewT Sounding
def skewT (t_temp,td_temp,tv_temp,wetzero,p_temp,u,v,parcel,preparcel,lclp,lclt,height,rstu, rstv, lstu, lstv, mw6, 
            cape, wmaxshr, q_mean, lapse, brn, name5,wenn,titel):

    register_projection(SkewXAxes)

    #fig = plt.figure(figsize=(6.5875, 6.2125))
    fig = plt.figure(figsize=(9, 6.2125))
    ax = fig.add_subplot(111, projection='skewx')
    plt.subplots_adjust(left=0.11, right=0.6, bottom=0.07, top=0.91)
    #Grid
    ax.grid(True)
    #lines
    pmax = 1000
    pmin = 100
    dp = -10
    presvals = np.arange(pmax, pmin+dp, dp)

    # plot the moist-adiabats
    for t in np.array([-40, -30, -20, -10,  0,  5, 10, 15, 20, 25, 30, 35, 40]):
        tw = []
        for p in presvals:
            tw.append(wetlift(1000., t, p))
        ax.semilogy(tw, presvals, 'b--', lw=0.5, alpha=0.6)

    # plot the dry adiabats
    for t in np.array([-80, -60, -40, -30, -20, -10,   0,  10,  20,  30,  40, 50,  70,  90, 110, 130, 150, 170, 190]):
        ax.semilogy(thetas(t, presvals), presvals, 'r--',lw=0.5, alpha=0.6)

    # linie gleichensaettiungsmischungsverhaeltnis
    pmin = 700
    presvals = np.arange(pmax, pmin+dp, dp)

    mrr = np.array([0.1,1.0,5.0,10.0,20.0,30.0])
    for i in range(0,np.size(mrr)):
        ax.semilogy(temp_at_mixrat(mrr[i], presvals),presvals,'m--',lw=0.6, alpha=0.7)

    #plot tile/location
    #titel = "Forecast Sounding at " + str(lat) + "N  " + str(lon) + "E"
    plt.annotate(titel, xy=(0.25, 1.04), xycoords='axes fraction',fontsize=15)
    plt.annotate(wenn, xy=(0, 1.005), xycoords='axes fraction',fontsize=10)
    plt.xlabel("Temperature [C]")
    plt.ylabel("Pressure [hPa]")

    l = ax.axvline(0, color='b', linestyle='--' , lw=0.9)
    l = ax.axvline(-20, color='b', linestyle='--', lw=0.9)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    ax.semilogy(wetzero, p_temp, 'c-', lw=0.8) # Plot the wetbulb profile
    ax.semilogy(tv_temp,p_temp,'b-'  , lw=1.1, alpha=0.9) # Plot the virtuell temperature
    ax.semilogy(t_temp, p_temp, 'b-' , lw=1.8 ) # Plot the temperature profile
    ax.semilogy(td_temp, p_temp, 'g-', lw=1.8) # plot the dewpoint profile
    ax.semilogy(parcel,preparcel,'r--',lw=1.4)

    mu = np.zeros(p_temp.size, dtype=np.float64)
    for i in range(0,p_temp.size):
        #if p_temp[i] >= lclp:
        #    mu[i] = thetas(theta(p_temp[0], tv_temp[0], p2=1000.), p_temp[i])
        #else:
        if p_temp[i] <= lclp:
            mu[i] = wetlift(lclp,lclt,p_temp[i])
    ax.fill_betweenx(p_temp,tv_temp, mu,  where=mu > tv_temp, facecolor='red', interpolate=True, alpha=0.4)

    # Disables the log-formatting that comes with semilogy
    ax.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_yticks(np.linspace(100,1000,10))
    ax.set_ylim(1050,100)
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.set_xlim(-30,40)

    #plt.annotate('-10', xy=(0.3, -0.04), xycoords='axes fraction',fontsize=10)
    #plt.annotate('-20', xy=(0.2, -0.04), xycoords='axes fraction',fontsize=10)
    #plt.annotate('-30', xy=(0.1, -0.04), xycoords='axes fraction',fontsize=10)

    # Draw the hodograph on the Skew-T.
    """
    ax2 = plt.axes([.685,.655,.25,.25])
    i=0
    #1km
    m=[]
    n=[]
    while height[i] < 1000 :
	    m.append(u[i])
	    n.append(v[i])
	    i += 1

    m.append(u[i])
    n.append(v[i])
    ax2.plot(m, n, 'g-', lw=1.5, label='0-1 km')

    m=[]
    n=[]
    while height[i] < 8000 :
            m.append(u[i])
            n.append(v[i])
            i += 1

    ax2.plot(m, n, 'k-', lw=1.5)
    """
    ax2 = plt.axes([.625,.40,.25,.50])

    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    clim=20
    if clim > 20:
        for i in range(10,clim+20,10):
            # Draw the range rings around the hodograph.
            circle = plt.Circle((0,0),i,color='k',alpha=.3, fill=False)
            ax2.add_artist(circle)
    else:
        for i in range(5,40,5):
            # Draw the range rings around the hodograph.
            circle = plt.Circle((0,0),i,color='k',alpha=.3, fill=False)
            ax2.add_artist(circle)

    i=0
    while height[i] < 10000 :
            i += 1
    points = np.array([u[:i+1], v[:i+1]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # Create a continuous norm to map from data points to colors
    #cmap = ListedColormap(['r', 'g', 'b'])
    norm = plt.Normalize(0, 10000)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","purple","navy","forestgreen","darkgreen"])
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    # Set the values used for colormapping
    lc.set_array(height)
    lc.set_linewidth(2)
    line = ax2.add_collection(lc)
    cbar=fig.colorbar(line, ax=ax2, orientation='horizontal', ticks=np.arange(0,10000,1000))
    clevs=np.arange(0,10,1)
    cbar.ax.set_xticklabels(labels=clevs[::1])

    ax2.plot(rstu, rstv, 'ro',alpha=0.9,zorder=2) # Plot Bunker's Storm motion right mover as a red dot
    ax2.plot(lstu, lstv, 'bo',alpha=0.9,zorder=2) # Plot Bunker's Storm motion left mover as a blue dot
    ax2.plot(mw6[0], mw6[1], 'o', color = 'grey',alpha=0.7,zorder=2) # 0-6km Mean weight wind

    ax2.axhline(y=0, color='k',zorder=1,lw=1.1)
    ax2.axvline(x=0, color='k',zorder=1,lw=1.1) 

    ax2.set_xlim(-clim,clim)
    ax2.set_ylim(-clim,clim)
    if clim == 40:
        ax2.text(-5.2,-39.7,"40 m/s")
        ax2.text(-5.2,-29.7,"30")
        ax2.text(-5.2,-19.7,"20")
        #ax2.text(-4.82,-9.7,"10 m/s")
    elif clim == 20:
        #ax2.text(-4.43,-29.7,"30 m/s")
        ax2.text(-3.6,-19.7,"20 m/s")
        ax2.text(-3.6,-9.7,"10 m/s")
    elif clim == 50:
        ax2.text(-6.2,-49.7,"50 m/s")
        ax2.text(-6.2,-39.7,"40")
        ax2.text(-6.2,-29.7,"30")
    else:
        ax2.text(-4.48,-29.7,"30 m/s")
        ax2.text(-4.48,-19.7,"20 m/s")
        ax2.text(-4.48,-9.7,"10 m/s")
    #ax2.annotate("10 m/s", xy=(0.5, 0.5), xycoords='axes fraction',fontsize=15)
    #fig.text(0.750, 0.38, "km")
    fig.text(0.750, 0.38, "km")
    
    val = 180./np.pi
    direction = np.arctan2(-mw6[0],-mw6[1]) * val
    if ( direction < 0):
        direction += 360;
    fig.text(0.625, 0.33, r"Mean Wind (0-6km) : %3.0f$^{\circ}$/%.0f m/s" % (direction,np.sqrt(mw6[0]*mw6[0]+mw6[1]*mw6[1])))
    """
    direction = np.arctan2(-mw8[0],-mw8[1]) * val
    if ( direction < 0):
        direction += 360;
    fig.text(0.625, 0.27, r"Mean Wind (0-8km) : %3.0f$^{\circ}$/%.0f m/s" % (direction,np.sqrt(mw8[0]*mw8[0]+mw8[1]*mw8[1])))
    """
    #fig.text(0.625, 0.30, "ML CAPE : %.1f J/kg" % (2843.2))
    fig.text(0.625, 0.27, "MU CAPE : %.1f J/kg" % (cape))
    fig.text(0.625, 0.24, "BRN ML : %.1f" % (brn))
    if (wmaxshr >= 1000.0):
        fig.text(0.625,0.21, r"WMAXSHEAR ML : %.1f $m^2/s^2$" % (wmaxshr),color='r',weight='bold')
    elif (wmaxshr < 1000.0) and (wmaxshr >= 400.0):
        fig.text(0.625,0.21, r"WMAXSHEAR ML : %.1f $m^2/s^2$" % (wmaxshr),color='orange',weight='bold')
    else:
        fig.text(0.625,0.21, r"WMAXSHEAR ML : %.1f $m^2/s^2$" % (wmaxshr)) 
    """
    if pw >= 40.0:
        fig.text(0.625,0.15,r"PW : %.1f $kg/m^2$" % (pw),color='r',weight='bold')
    elif (pw < 40.0) and (pw > 30.0):
        fig.text(0.625,0.15,r"PW : %.1f $kg/m^2$" % (pw),color='orange',weight='bold')
    else:
        fig.text(0.625,0.15,r"PW : %.1f $kg/m^2$" % (pw)) 
    """
    fig.text(0.625, 0.18, "Specific humidity (100mb) : %.1f g/kg" % (q_mean))
    fig.text(0.625, 0.15, "Lapse Rate (850-500hPa) : %.1f K/km" %(lapse))

    try:
        plt.savefig(name5)
    except ValueError:
        pass
    plt.close()

skewT (t_temp,td_temp,vir_temp,wetzero,p_temp,u,v,mu_parcel,press_parcel,float(sys.argv[4]),
       float(sys.argv[5]),height,float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]), float(sys.argv[9]), (mwu6, mwv6),
       float(sys.argv[3]),float(sys.argv[12]), float(sys.argv[15]), float(sys.argv[14]), float(sys.argv[13]), name5,wenn,title)
