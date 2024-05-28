# ---------------------------------------------------------------------------------------------------------------------
#                  Entraining CAPE Library
#
# https://arxiv.org/pdf/2301.04712
# https://figshare.com/articles/software/ECAPE_scripts/21859818?file=42303630
#
#
# ---------------------------------------------------------------------------------------------------------------------


import numpy as np
from scipy.special import lambertw
from scipy.special import sici

# ---------------------------------------------------------------------------------------------------------------------

def comp_cdwave(F):
   gamma_em = np.euler_gamma
   return 4*(F*np.sin(2/F) - (F**2)*(np.sin(1/F)**2) - cosint(2/F) + np.log(2/F) - 1 + gamma_em)


def cosint(x):
   _, ci = sici(x)
   return ci

# ---------------------------------------------------------------------------------------------------------------------
# descriminator function between liquid and ice (i.e., omega defined in the
# beginning of section 2e in Peters et al. 2022)

def omega(T,T1,T2):
    return ((T - T1)/(T2-T1))*np.heaviside((T - T1)/(T2-T1), 1)*np.heaviside((1 - (T - T1)/(T2-T1)), 1) + np.heaviside(-(1 - (T - T1)/(T2-T1)), 1)


def domega(T,T1,T2):
    return (np.heaviside(T1-T, 1) - np.heaviside(T2-T, 1))/(T2-T1)

# ---------------------------------------------------------------------------------------------------------------------
