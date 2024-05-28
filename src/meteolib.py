#!/usr/bin/python3
import math
import numpy as np

# ---------------------------------------------------------------------------------------------------------------------
# Source for the most functions are sharppy

# ---------------------------------------------------------------------------------------------------------------------
# meteorology constants

cr = {'ZEROCNK': 273.15,                       # K        Zero Celsius in Kelvins
      'Rd': 287.04,                            # J/kg/K   dry air gas constant
      'Rv': 461.5,                             # J/kg/K   water vapor gas constant
      'cpd': 1005.7,                           # J/kg/K   isobaric dry air specific heat capacit
      'cpv': 1850.0,   # 1870                  # J/kg/K   water vapor specific heat capacity (isobaric),
      'cl': 4218.0,    # cpl 4190              # J/kg/K   liquid water specific heat capacity
      'cpi': 2106.0,                           # J/kg/K   liquid water specific heat capacity
      'eps': 0.62197,                          # dimless  molecular weight ratio
      'ROCP': 0.28571426,                      # dimless  R over Cp
      'Lref': 2.5*1e6,                         # J/kg     latent heat of condensation at reference conditions
      'pref': 611.65,                          # Pa       partial pressure of water vapor a triple point temp
      'ttrip': 273.15,                         # K        triple point temp
      'G': 9.80665,                            # m/(s*s)  Gravity
      'Re': 6371008.767,                       # m        earth radius
      # 'omega': (2*math.pi)/(23*3600+56*60+4.1),   # rad/s    angular velocity of the earth
      'omega': 7.29211505392569e-05,           # rad/s    angular velocity of the earth
      'p0': 100000.0,                          # Pa       reference pressure
      'T0': 273.15,                            # K        reference temperature
      'e0': 611,                               # Pa       reference water vapor (p0, T0)
      'kn2ms': 0.514,                          # 1 knot = 0.514 m/s
      'ms2kn': 1.94384,                        # 1 m/s  = 1.94384
      'SeaMile2m': 1852.,                      # Sea Miles to m
      'StatuteMile2m': 1609.,                  # Statue Mile to m
      'Pi': math.pi,                           # Pi
      'e': 2.71828,                            # Euler Zahl
      'Pas2hPah': 36,                          # Pa/s to hPa/h
      'ksq': 0.18,                             # VON KARMAN CONSTANT
      'Pr': 1/3,                               # PRANDTL NUMBER
      'Lm': 250,                               # HORIZONTAL MIXING LENGTH
      'c1': 0.0498646455,
      'c2': 2.4082965,
      'c3': 7.07475,
      'c4': 38.9114,
      'c5': 0.0915,
      'c6': 1.2035
      }

# ---------------------------------------------------------------------------------------------------------------------


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
            pwrp = np.power((p / 1000.), cr['ROCP'])
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
    return 1000. / (np.power((theta / t), (1./cr['ROCP'])))


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
    return (temp * np.power((p0 / pres), cr['ROCP']))


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
    Cp/R ln(T0/T) = p0/p |T0 = Theata
    p = p0 / (T0/T)^Cp/R
    '''
    thalvl = p0 / (pow((theta / temp), (1/cr['ROCP'])))
    return thalvl


def thetas(theta, presvals, p0=1000.):
    '''
    Parameters
    ----------
    Returns
    -------
    Temperature in K as numpy.array
    '''
    return (theta / (np.power((p0 / presvals), cr['ROCP'])))


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

# Source: SHARPpy
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
    return (t + ((2.501*1000000)/1004)*q) * np.power((p0/p), cr['ROCP'])

# ---------------------------------------------------------------------------------------------------------------------


def uvwind(winddir, wind_speed):
    """
    This function convert wind speed and direction in to u,v component.
    The Winddirection is given in degree,
    """

    u = (-1) * np.multiply(wind_speed, np.sin(np.radians(winddir)))
    v = (-1) * np.multiply(wind_speed, np.cos(np.radians(winddir)))
    return u, v


def uv2spddir(u, v):

    direction = np.rad2deg(np.arctan2(-u, -v))

    if isinstance(direction, np.ndarray):
        direction = np.remainder(direction + 360, 360)
    else:
        direction = (direction + 360) % 360

    return (direction, np.sqrt(u*u + v*v))


def mean_wind(u, v, ps, stu=0.0, stv=0.0):
    return np.average(u, weights=ps) - stu, np.average(v, weights=ps) - stv


def shear(ubot, vbot, utop, vtop):
    du = utop - ubot
    dv = vtop - vbot
    return np.sqrt(du*du + dv*dv)
