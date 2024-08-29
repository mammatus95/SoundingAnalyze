# ----------------------------------------------------------------------------------------------------------------------------
#                  Entraining CAPE Library
#
# https://arxiv.org/pdf/2301.04712
# https://figshare.com/articles/software/ECAPE_scripts/21859818?file=42303630
#
#
# ----------------------------------------------------------------------------------------------------------------------------
# NOTE: most of the scripts and functions that use this function need
# saturation mass fraction qs, not saturation mixing ratio rs.  To get
# qs from rs, use the formula qs = (1 - qt)*rs, where qt is the total
# water mass fraction

import numpy as np
from scipy.special import lambertw
from scipy.special import sici
from meteolib import cr
# ----------------------------------------------------------------------------------------------------------------------------

def comp_cdwave(F):
   gamma_em = np.euler_gamma
   return 4*(F*np.sin(2/F) - (F**2)*(np.sin(1/F)**2) - cosint(2/F) + np.log(2/F) - 1 + gamma_em)


def cosint(x):
   _, ci = sici(x)
   return ci

# ----------------------------------------------------------------------------------------------------------------------------
# descriminator function between liquid and ice (i.e., omega defined in the
# beginning of section 2e in Peters et al. 2022)

def omega(T,T1,T2):
    return ((T - T1)/(T2-T1))*np.heaviside((T - T1)/(T2-T1), 1)*np.heaviside((1 - (T - T1)/(T2-T1)), 1) + np.heaviside(-(1 - (T - T1)/(T2-T1)), 1)


def domega(T,T1,T2):
    return (np.heaviside(T1-T, 1) - np.heaviside(T2-T, 1))/(T2-T1)

# ----------------------------------------------------------------------------------------------------------------------------
def get_qs(qt, rs):
    """
    Calculate the saturation mass fraction qs from the total water mass fraction qt and the saturation mixing ratio rs.

    Parameters
    ----------
    qt : numpy.ndarray
        Total water mass fraction (kg/kg)
    rs : numpy.ndarray
        Saturation mixing ratio (kg/kg)

    Returns
    -------
    qs : numpy.ndarray
        Saturation mass fraction (kg/kg)
    """
    return (1 - qt)*rs

# ----------------------------------------------------------------------------------------------------------------------------
# saturation mixing ratio
def compute_rsat(T, p, T1, T2, iceflag=0):

    omeg = omega(T,T1,T2)
    epsilon = cr['Rd']/cr['Rv']
    
    if iceflag == 0:
        term1 = (cr['cpv']-cr['cpl'])/cr['Rv']
        term2 = (cr['xlv']-cr['ttrip']*(cr['cpv']-cr['cpl']))/cr['Rv']
        esl = np.exp((T-cr['ttrip'])*term2/(T*cr['ttrip']))*cr['eref']*(T/cr['ttrip'])**(term1)
        qsat = epsilon*esl/(p-esl)
    elif iceflag == 1: # give linear combination of mixing ratio with respect to liquid and ice (eq. 20 in Peters et al. 2022)
        term1 = (cr['cpv']-cr['cpl'])/cr['Rv']
        term2 = (cr['xlv']-cr['ttrip']*(cr['cpv']-cr['cpl']))/cr['Rv']
        esl_l = np.exp((T-cr['ttrip'])*term2/(T*cr['ttrip']))*cr['eref']*(T/cr['ttrip'])**(term1)
        qsat_l = epsilon*esl_l/(p-esl_l)
        term1 = (cr['cpv']-cr['cpl'])/cr['Rv']
        term2 = (cr['xls']-cr['ttrip']*(cr['cpv']-cr['cpl']))/cr['Rv']
        esl_i = np.exp((T-cr['ttrip'])*term2/(T*cr['ttrip']))*cr['eref']*(T/cr['ttrip'])**(term1)
        qsat_i = epsilon*esl_i/(p-esl_i)
        qsat = (1-omeg)*qsat_l + (omeg)*qsat_i
    elif iceflag == 2: # only give mixing ratio with respect to ice
        term1 = (cr['cpv']-cr['cpi'])/cr['Rv']
        term2 = (cr['xls']-cr['ttrip']*(cr['cpv']-cr['cpi']))/cr['Rv']
        esl = np.exp((T-cr['ttrip'])*term2/(T*cr['ttrip']))*cr['eref']*(T/cr['ttrip'])**(term1)
        esl = min( esl , p*0.5 )
        qsat = epsilon*esl/(p-esl)
    return qsat


# ----------------------------------------------------------------------------------------------------------------------------
# lapse rate for an unsaturated parcel
def drylift(T, qv, T0, qv0, fracent):
    
    """
    Calculate the lapse rate of an unsaturated parcel.

    Parameters
    ----------
    T : float
        Temperature of parcel in K
    qv : float
        Specific humidity of parcel in kg/kg
    T0 : float
        Temperature at the reference level in K
    qv0 : float
        Specific humidity at the reference level in kg/kg
    fracent : float
        fractional entrainment rate (in m^-1)

    Returns
    -------
    gamma_d : float
        lapse rate of the parcel in K/m
    """
    cpmv = (1 - qv)*cr['cpd'] + qv*cr['cpv']
    B = cr['G']*( (T-T0)/T0 + (cr['Rv']/cr['Rd'] - 1)*(qv - qv0) )
    eps = -fracent*(T - T0)
    gamma_d = - (cr['G'] + B)/cpmv + eps
    return gamma_d

# ----------------------------------------------------------------------------------------------------------------------------
# lifted condensation level
# using the romps 2017 formula

def compute_LCL(T, qv, p):
    cpm = (1 - qv)*cr['cpd'] + qv*cr['cpv']
    Rm = (1 - qv)*cr['Rd'] + qv*cr['Rv']
    
    a = cpm/Rm + ( cr['cl'] - cr['cpv'] )/cr['Rv']
    b = -(cr['xlv'] - (cr['cpv'] - cr['cl'])*cr['ttrip'])/(cr['Rv']*T)
    c = b/a
    
    r_sat = compute_rsat(T,p,0,273.15,253.15)
    q_sat = r_sat/(1 + r_sat)
    RH = qv/q_sat
    arg1 = RH**(1/a)
    arg2 = c*np.exp(1)**c
    arg3 = lambertw(arg1*arg2,k=-1)
    T_LCL = c*T/arg3
    P_LCL = p*(T_LCL/T)**(cpm/Rm)
    Z_LCL = (cpm/cr['G'])*(T - T_LCL)
    
    return Z_LCL


# using numerical integration
# using numerical integration
def compute_LCL_NUMERICAL(T,qv,p,dz):

    nfound_LCL = True
    ind_hgt = 0
    Ton = T
    Qon = qv
    Pon = p
    while nfound_LCL:
        ind_hgt = ind_hgt+1
        Ton = Ton + dz*drylift(Ton,Qon,Ton,Qon,0)
        Pon = Pon - dz*(Pon*cr['G'])/(cr['Rd']*(1 + (cr['Rv']/cr['Rd'] - 1)*Qon )*Ton )
        rsat = compute_rsat(Ton,Pon,0,273.15,253.15)
        qsat = rsat/(1 + rsat)
        if Qon >= qsat:
            nfound_LCL = False
    Z_LCL = ind_hgt*dz
    
    return Z_LCL

# ----------------------------------------------------------------------------------------------------------------------------
# lapse rate for a saturated parcel

def moislif(T, qv, qvv, qvi, p0, T0, q0, qt, fracent, prate, T1, T2):

    epsilon = cr['Rd']/cr['Rv']
 
    qt = max(qt,0.0)
    qv = max(qv,0.0)
    
    OMEGA = omega(T,T1,T2)
    dOMEGA = domega(T,T1,T2)
    
    
    cpm = (1 - qt)*cr['cpd'] + qv*cr['cpv'] + (1 - OMEGA)*(qt-qv)*cr['cl'] + OMEGA*(qt-qv)*cr['cpi']
    Lv = cr['xlv'] + (T - cr['ttrip'])*(cr['cpv'] - cr['cl'])
    Li = (cr['xls']-cr['xlv']) + (T - cr['ttrip'])*(cr['cl'] - cr['cpi'])
    Rm0 = (1 - q0)*cr['Rd'] + q0*cr['Rv']
    

    T_rho=T*(1 - qt + qv/epsilon)
    T_rho0=T0*( 1 - q0 + q0/epsilon )
    B = g*(T_rho - T_rho0)/(T_rho0)
    
    Qvsl = qvv/( epsilon - epsilon*qt + qv)
    Qvsi = qvi/( epsilon - epsilon*qt + qv)
    Q_M = (1 - OMEGA)*qvv/(1 - Qvsl) + OMEGA*qvi/(1 - Qvsi)
    L_M = Lv*(1 - OMEGA)*qvv/(1 - Qvsl) + (Lv + Li)*OMEGA*qvi/(1 - Qvsi)

    
    eps_T = -fracent*(T - T0)
    eps_qv = -fracent*(qv - q0)
    eps_qt = -fracent*(qt - q0)-prate*(qt-qv)
    term1 = -B
    
    term2 = - Q_M*(Lv + Li*OMEGA)*g/(Rm0*T0)
    
    term3 = -g
    term4 = (cpm - Li*(qt - qv)*dOMEGA)*eps_T
    term5 = (Lv + Li*OMEGA)*(eps_qv + (qv/(1-qt))*eps_qt)

    term6 = cpm
    term7 = -Li*(qt - qv)*dOMEGA
    term8 = (Lv + Li*OMEGA)*(-dOMEGA*(qvv - qvi) + (1/(cr['Rv']*(T**2)))*(L_M))
    gamma_m = ( term1 + term2 + term3 + term4 + term5)/(term6 + term7 + term8)
    return gamma_m