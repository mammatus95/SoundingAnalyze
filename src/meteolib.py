#!/usr/bin/python3
import math
import numpy as np

import src.meteo_constants as meteo_constants
# -----------------------------------------------------------------------------------------------------------------------------
# Source for the most functions are sharppy

# -----------------------------------------------------------------------------------------------------------------------------
# meteorology constants

cr = meteo_constants.constants

# -----------------------------------------------------------------------------------------------------------------------------


def virtuelle(q, t):
    """
    Parameters:
    ----------
    w : specific humidity (g/kg) as numpy.array
    t : temperature (K) as numpy.array
    Returns:
    ----------
    virtual temperature (K) as numpy.array
    """
    return t*(1+0.608*q)


def q_to_mixrat(q):
    """
    Parameters
    ----------
    q : Specific Humidity (kg/kg) as numpy.array
    Returns:
    ----------
    Mixing Ratio (kg/kg) as numpy.array
    """
    return q/(1.0-q)


def mixrat_to_q(mr):
    """
    Parameters
    ----------
    q : Mixing Ratio (kg/kg) as numpy.array
    """
    return mr/(1.0+mr)


def temp_at_mixrat(w, p):
    '''
    Returns the temperature (C) of air at the given mixing ratio (g/kg) and
    pressure (hPa)
    Parameters
    ----------
    w : Mixing Ratio (g/kg) as numpy.array
    p : Pressure (hPa) as numpy.array
    Returns
    -------
    Temperature (C) of air at given mixing ratio and pressure
    '''

    x = np.log10(w * p / (622. + w))
    x = (np.power(10., ((cr['c1'] * x) + cr['c2'])) - cr['c3'] +
         (cr['c4'] * np.power((np.power(10, (cr['c5'] * x)) - cr['c6']), 2))) - cr['ZEROCNK']
    return x


def lcltemp(t, td):
    '''
    Parameters
    ----------
    t : Temperature of the parcel in C as float
    td : Dewpoint temperature of the parcel in C as float
    Returns
    -------
    Temperature of lcl level in C as float
    '''
    s = t - td
    dlt = s * (1.2185 + 0.001278 * t + s * (-0.00219 + 1.173e-5 * s -
        0.0000052 * t))
    return t - dlt


def lclhigh (t, td, high):
    return high + ((t-td) * 125.0)


def satdampfdruck (t):
    # wahrscheindlich der Saettigungsdampfdruck fur wasser also der maximale partialdruck den der wasserdampf erreichen kann
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
        npol = 15.13 / (np.power(npol, 4))
        ppol = t * (4.9618922e-07 + t * (-6.1059365e-09 + t * (3.9401551e-11 + t * (-1.2588129e-13 + t * (1.6688280e-16)))))
        ppol = 1 + t * (3.6182989e-03 + t * (-1.3603273e-05 + ppol))
        ppol = (29.93 / np.power(ppol, 4)) + (0.96 * t) - 14.8
        correction = np.zeros(t.shape, dtype=np.float64)
        correction[t <= 0] = npol[t <= 0]
        correction[t > 0] = ppol[t > 0]
        return correction
    else:
        if t is np.ma.masked:
            return t
        if t <= 0:
            npol = 1. + t * (-8.841660499999999e-3 + t * ( 1.4714143e-4 + t * (-9.671989000000001e-7 + t * (-3.2607217e-8 + t * (-3.8598073e-10)))))
            npol = 15.13 / (np.power(npol, 4))
            return npol
        else:
            ppol = t * (4.9618922e-07 + t * (-6.1059365e-09 + t * (3.9401551e-11 + t * (-1.2588129e-13 + t * (1.6688280e-16)))))
            ppol = 1 + t * (3.6182989e-03 + t * (-1.3603273e-05 + ppol))
            ppol = (29.93 / np.power(ppol, 4)) + (0.96 * t) - 14.8
            return ppol


def satlift(p, thetam):
    if np.fabs(p - 1000.) - 0.001 <= 0: return thetam
    eor = 999
    e2, t2 = np.nan, np.nan
    while np.fabs(eor) - 0.1 > 0:
        if eor == 999:                  # First Pass
            pwrp = np.power((p / 1000.), cr['kappa'])
            t1 = (thetam + cr['ZEROCNK']) * pwrp - cr['ZEROCNK']
            e1 = wobf(t1) - wobf(thetam)
            rate = 1
        else:                           # Successive Passes
            rate = (t2 - t1) / (e2 - e1)
            t1 = t2
            e1 = e2
        t2 = t1 - (e1 * rate)
        e2 = (t2 + cr['ZEROCNK']) / pwrp - cr['ZEROCNK']
        e2 += wobf(t2) - wobf(e2) - thetam
        eor = e2 * rate
    return t2 - eor


def thalvl(theta, t):
    '''
    Returns the level (hPa) of a parcel.
    Parameters
    ----------
    theta : Potential temperature of the parcel in K as numpy.array
    t : Temperature of the parcel in K as numpy.array

    Returns
    -------
    Pressure level in hPa
    '''
    return 1000. / (np.power((theta / t), (1./cr['kappa'])))


def theta(pres, temp, p0=1000.):
    '''
    Parameters
    ----------
    pres : The pressure of the parcel in hPa as numpy.array
    temp : Temperature of the parcel in K as numpy.array
    p0   : Reference pressure level in hPa as float (default 1000.)
    Returns
    -------
    Potential temperature in K
    '''
    return (temp * np.power((p0 / pres), cr['kappa']))


def potlvl(theta, temp, p0=1000.):
    '''
    Parameters
    ----------
    theta : Potential temperature of the parcel in K as numpy.array
    temp : Temperature of the parcel in K as numpy.array

    Returns
    -------
    Pressure Level (hPa [float]) of the parcel

    Cp (ln T - ln T0) = R (ln p - ln p0)
    Cp/R ln(T/T0) = ln(p/p0) |T0 = Theata
    p/p0 = (T/T0)^Cp/R
    p = p0 / (T/T0)^Cp/R
    '''
    thalvl = p0 / (pow((theta / temp), (1./cr['kappa'])))
    return thalvl


def thetas(theta, presvals, p0=1000.):
    '''
    Parameters
    ----------
    Returns
    -------
    Temperature in K as numpy.array

    Cp (ln T - ln T0) = R (ln p - ln p0)
    Cp/R ln(T/T0) = ln(p/p0)
    T = T0 * (p/p0)^R/Cp    
    '''
    return (theta * np.power(presvals/p0, cr['kappa']))


def drylift(p, t, td):
    '''
    Lifts a parcel to the LCL and returns its new level and temperature.
    Parameters
    ----------
    p : Pressure of initial parcel in hPa as float
    t : Temperature of inital parcel in K as float
    td : Dew Point of initial parcel in K as float
    Returns
    -------
    p2 : LCL pressure as float
    t2 : LCL Temperature as float
    '''
    t2 = lcltemp(t - cr['ZEROCNK'], td - cr['ZEROCNK']) + cr['ZEROCNK']
    p2 = thalvl(theta(p, t, 1000.), t2)
    return p2, t2


def wetlift(p, t, p2):
    '''
    Lifts a parcel to the LCL and returns its new level and temperature.
    Parameters
    ----------
    p : Pressure of initial parcel in hPa as float
    t : Temperature of inital parcel in K as float
    p2:
    Returns
    -------
    Temperature in K as numpy.array
    '''
    thta = theta(p, t, 1000.) - cr['ZEROCNK']
    thetam = thta - wobf(thta) + wobf(t - cr['ZEROCNK'])
    return satlift(p2, thetam) + cr['ZEROCNK']


def wetpot(p, t, p2=1000.):
    return wetlift(p, t, p2)


def wetbulb(p, t, td):
    '''
    Parameters
    ----------
    p : Pressure of parcel in hPa as float
    t : Temperature of parcel in K as float
    td : Dew Point of parcel in K as float
    Returns
    -------
    Wetbulb temperature in K as float
    '''
    p2, t2 = drylift(p, t - cr['ZEROCNK'], td - cr['ZEROCNK'])
    return wetlift(p2, t2, p) + cr['ZEROCNK']


def thetae(p, t, q, p0=1000.):
    """
    Calculates the equlivalent potential temperature (K) for the given parcel
    Parameters
    ----------
    p : Pressure of parcel in hPa as numpy.array
    t : Temperature of parcel in K as numpy.array
    q : specific humidity in kg/kg as numpy.array
    Returns
    -------
    qulivalent potential temperature in K as numpy.array
    """
    return (t + ((cr['Lref'])/cr['cpd'])*q) * np.power((p0/p), cr['kappa'])

# -----------------------------------------------------------------------------------------------------------------------------


def uvwind(winddir, wind_speed):
    """
    This function convert wind speed and direction in to u,v component.
    The Winddirection is given in degree,
    """
    winddir, wind_speed = np.abs(winddir), np.abs(wind_speed)
    u = (-1) * np.multiply(wind_speed, np.sin(np.radians(winddir)))
    v = (-1) * np.multiply(wind_speed, np.cos(np.radians(winddir)))
    return u, v


def uv2spddir(u, v):
    direction = np.rad2deg(np.arctan2(-u, -v))
    if isinstance(direction, np.ndarray):
        direction = np.remainder(direction + 360, 360)
    else:
        direction = (direction + 360) % 360

    wind_speed = np.sqrt(np.square(u) + np.square(v))
    if type(wind_speed) is not np.ndarray:
        if wind_speed == 0:
            direction = np.nan
    else:
        direction[np.where(wind_speed == 0)] = np.nan

    return (direction, wind_speed)


def mean_wind(u, v, ps, stu=0.0, stv=0.0):
    return np.average(u, weights=ps) - stu, np.average(v, weights=ps) - stv


def shear(ubot, vbot, utop, vtop):
    du = utop - ubot
    dv = vtop - vbot
    return np.sqrt(du*du + dv*dv)
