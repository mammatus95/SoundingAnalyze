# ----------------------------------------------------------------------------------------------------------------------------
#                  Entraining CAPE Library
# Author : John M. Petersa
# Rework : Morten Kretschmer
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
from src.meteolib import cr
#from wavefunction_lib import comp_cdwave



# ----------------------------------------------------------------------------------------------------------------------------
# descriminator function between liquid and ice (i.e., omega defined in the
# beginning of section 2e in Peters et al. 2022)

def omega(T, T1, T2):
    normalized_T = (T - T1) / (T2 - T1)
    return np.clip(normalized_T, 0, 1)

def domega(T, T1, T2):
    """
    Compute the derivative of the discriminator function omega(T, T1, T2).

    - When T < T1, both heaviside functions are 1, so the result is 0.
    - When T1 <= T < T2, the first heaviside function (T1 - T) is 0 and the second (T2 - T) is 1, so the result is -1.
    - When T >= T2, both heaviside functions are 0, so the result is 0.

    Parameters
    ----------
    T : numpy.ndarray
        Temperature in K
    T1 : float
        Warmest mixed-phase temperature
    T2 : float
        Coldest mixed-phase temperature

    Returns
    -------
    d_omega : numpy.ndarray
        Derivative of omega with respect to T

    Notes
    -----
    The derivative is computed using the Heaviside step function, which is
    zero for T < T1 and T > T2, and one for T1 <= T <= T2.  The derivative
    is then computed as the difference of the Heaviside functions divided by
    T2-T1.  This is a simple approximation, which should be sufficient for
    most purposes.  If a more accurate derivative is needed, it should be
    computed numerically using a more sophisticated method.
    """
    # T1 equal to T2 should raise a ZeroDivisionError due to rounding errors (probably not intended by author)

    if T1 == T2:
        raise ZeroDivisionError("T1 is equal T2")
    
    mask = (T2 <= T) & (T <= T1)
    return np.where(mask, 1/(T2 - T1), 0)

# ----------------------------------------------------------------------------------------------------------------------------


def get_qs(qt, rs):
    """
    Calculate the saturation mass fraction qs from the total water mass fraction qt and the saturation mixing ratio rs.

    Parameters
    ----------
    qt : numpy.ndarray
        Total water mass fraction in kg/kg
    rs : numpy.ndarray
        Saturation mixing ratio in kg/kg

    Returns
    -------
    qs : numpy.ndarray
        Saturation mass fraction in kg/kg
    """
    return (1 - qt)*rs


# ----------------------------------------------------------------------------------------------------------------------------
# moist static energy

def compute_moist_static_energy(T, q, z):
    """
    Compute the moist static energy (MSE) of a parcel.

    Parameters
    ----------
    T : numpy.ndarray
        Temperature of the parcel in K
    q : numpy.ndarray
        Total water mass fraction of the parcel in kg/kg
    z : numpy.ndarray
        Height above ground level of the parcel in m

    Returns
    -------
    MSE : numpy.ndarray
        Moist static energy of the parcel in J/kg

    Units
    -----
    [MSE] = J/kg/K * K + J/kg* kg/kg + m/(s*s) * m = [J/kg]
    """
    return cr['cpl']*T + cr['xlv']*q + cr['G']*z

# ----------------------------------------------------------------------------------------------------------------------------
# saturation mixing ratio


def compute_rsat(T, p, T1, T2, iceflag=0):
    """
    This function computes the saturation mixing ratio, using the integrated
    CLAUSIUS CLAPEYRON Equation (eq. 7-12 in Peters et al. 2022).
    
    Parameters
    ----------
    T : numpy.ndarray
        Temperature of the air (K)
    p : numpy.ndarray
        Pressure of the air (Pa)
    T1 : float
        warmest mixed-phase temperature?
    T2 : float
        coldest mixed-phase temperature?
    iceflag : int
        0: compute saturation mixing ratio with respect to liquid water
        1: compute saturation mixing ratio as a linear combination of mixing
            ratio with respect to liquid and ice water
        2: compute saturation mixing ratio with respect to ice

    Returns
    -------
    qsat : numpy.ndarray
        Saturation mixing ratio of the air (kg/kg)
    
    References
    ----------
    .. [1]  https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1

    """
    
    epsilon = cr['Rd']/cr['Rv']

    if iceflag == 0:
        term1 = (cr['cpv']-cr['cpl'])/cr['Rv']
        term2 = (cr['xlv']-cr['ttrip']*(cr['cpv']-cr['cpl']))/cr['Rv']
        esl_l = np.exp((T-cr['ttrip'])*term2/(T*cr['ttrip']))*cr['eref']*(T/cr['ttrip'])**term1
        return epsilon*esl_l/(p-esl_l)
    elif iceflag == 1: # give linear combination of mixing ratio with respect to liquid and ice (eq. 20 in Peters et al. 2022)
        omeg = omega(T, T1, T2)

        term1 = (cr['cpv']-cr['cpl'])/cr['Rv']
        term2 = (cr['xlv']-cr['ttrip']*(cr['cpv']-cr['cpl']))/cr['Rv']
        esl_l = np.exp((T-cr['ttrip'])*term2/(T*cr['ttrip']))*cr['eref']*(T/cr['ttrip'])**term1
        qsat_l = epsilon*esl_l/(p-esl_l)
        term1 = (cr['cpv']-cr['cpi'])/cr['Rv']
        term2 = (cr['xls']-cr['ttrip']*(cr['cpv']-cr['cpi']))/cr['Rv']
        esl_i = np.exp((T-cr['ttrip'])*term2/(T*cr['ttrip']))*cr['eref']*(T/cr['ttrip'])**term1
        qsat_i = epsilon*esl_i/(p-esl_i)
        return (1-omeg)*qsat_l + (omeg)*qsat_i
    elif iceflag == 2: # only give mixing ratio with respect to ice
        term1 = (cr['cpv']-cr['cpi'])/cr['Rv']
        term2 = (cr['xls']-cr['ttrip']*(cr['cpv']-cr['cpi']))/cr['Rv']
        esl_i = np.exp((T-cr['ttrip'])*term2/(T*cr['ttrip']))*cr['eref']*(T/cr['ttrip'])**term1
        esl_i = min(esl_i, p*0.5)
        return epsilon*esl_i/(p-esl_i)
    else:
        raise ValueError("iceflag must be 0, 1, or 2!")

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
        sounding profile of water vapor mass fraction (specific humidity) in kg/kg
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
    """
    Compute the lifted condensation level (LCL) following Romps (2017)

    Parameters
    ----------
    T : float
        Temperature in K
    qv : float
        Water vapor mass fraction (specific humidity) in kg/kg
    p : float
        Pressure in Pa

    Returns
    -------
    Z_LCL : float
        LCL height in m
    T_LCL : float
        LCL temperature in K
    P_LCL : float
        LCL pressure in Pa
    """
    cpm = (1 - qv)*cr['cpd'] + qv*cr['cpv']
    Rm = (1 - qv)*cr['Rd'] + qv*cr['Rv']
    
    a = cpm/Rm + ( cr['cpl'] - cr['cpv'] )/cr['Rv']
    b = -(cr['xlv'] - (cr['cpv'] - cr['cpl'])*cr['ttrip'])/(cr['Rv']*T)
    c = b/a
    
    r_sat = compute_rsat(T, p, 273.15, 253.15, 0)
    q_sat = r_sat/(1 + r_sat)
    RH = qv/q_sat
    arg1 = RH**(1/a)
    arg2 = c*np.exp(1)**c
    arg3 = np.real(lambertw(arg1*arg2, k=-1))
    T_LCL = c*T/arg3
    P_LCL = p*(T_LCL/T)**(cpm/Rm)
    Z_LCL = (cpm/cr['G'])*(T - T_LCL)
    
    return float(Z_LCL), float(T_LCL), float(P_LCL)


# using numerical integration
def compute_LCL_NUMERICAL(T, qv, p, dz):

    """
    Compute the lifted condensation level (LCL) using numerical integration

    Parameters
    ----------
    T : float
        Temperature in K
    qv : float
        Water vapor mass fraction (specific humidity) in kg/kg
    p : float
        Pressure in Pa
    dz : float
        Vertical resolution of the integration in m

    Returns
    -------
    Z_LCL : float
        LCL height in m
    """
    nfound_LCL = True
    ind_hgt = 0
    Ton = T
    Qon = qv
    Pon = p
    while nfound_LCL:
        ind_hgt = ind_hgt+1
        Ton = Ton + dz*drylift(Ton, Qon, Ton, Qon, 0)
        Pon = Pon - dz*(Pon*cr['G'])/(cr['Rd']*(1 + (cr['Rv']/cr['Rd'] - 1)*Qon )*Ton )
        rsat = compute_rsat(Ton, Pon, 273.15, 253.15, 0)
        qsat = rsat/(1 + rsat)
        if Qon >= qsat:
            nfound_LCL = False
    Z_LCL = ind_hgt*dz
    
    return float(Z_LCL)

# ----------------------------------------------------------------------------------------------------------------------------
# lapse rate for a saturated parcel

def moislif(T, qv, qvv, qvi, p0, T0, q0, qt, fracent, prate, T1, T2):

    """
    Calculate the lapse rate of a saturated parcel (gamma_m) using the 
    formulation from Peters et al. (2022).

    Parameters
    ----------
    T : float
        Temperature of parcel in K
    qv : float
        Water vapor mass fraction of parcel in kg/kg
    qvv : float
        Water vapor mass fraction of parcel in kg/kg, with respect to liquid
    qvi : float
        Water vapor mass fraction of parcel in kg/kg, with respect to ice
    p0 : float
        Pressure at the reference level in Pa
    T0 : float
        Temperature of sounding profile in K
    q0 : float
        Water vapor mass fraction (specific humidity) of sounding profile in kg/kg
    qt : float
        Total water mass fraction of parcel in kg/kg
    fracent : float
        fractional entrainment rate (in m^-1)
    prate : float
        precipitation rate (in m^-1)
    T1 : float
        warmest mixed-phase temperature
    T2 : float
        coldest mixed-phase temperature

    Returns
    -------
    gamma_m : float
        lapse rate of the parcel in K/m
    """

    epsilon = cr['Rd']/cr['Rv']
 
    qt = max(qt, 0.0)
    qv = max(qv, 0.0)
    
    OMEGA = omega(T, T1, T2)
    dOMEGA = domega(T, T1, T2)
    
    
    cpm = (1 - qt)*cr['cpd'] + qv*cr['cpv'] + (1 - OMEGA)*(qt-qv)*cr['cpl'] + OMEGA*(qt-qv)*cr['cpi']
    Lv = cr['xlv'] + (T - cr['ttrip'])*(cr['cpv'] - cr['cpl'])
    Li = (cr['xls']-cr['xlv']) + (T - cr['ttrip'])*(cr['cpl'] - cr['cpi'])
    Rm0 = (1 - q0)*cr['Rd'] + q0*cr['Rv']
    

    T_rho = T*(1 - qt + qv/epsilon)
    T_rho0 = T0*( 1 - q0 + q0/epsilon )
    B = cr['G']*(T_rho - T_rho0)/(T_rho0)
    
    Qvsl = qvv/( epsilon - epsilon*qt + qv)
    Qvsi = qvi/( epsilon - epsilon*qt + qv)
    Q_M = (1 - OMEGA)*qvv/(1 - Qvsl) + OMEGA*qvi/(1 - Qvsi)
    L_M = Lv*(1 - OMEGA)*qvv/(1 - Qvsl) + (Lv + Li)*OMEGA*qvi/(1 - Qvsi)

    
    eps_T = -fracent*(T - T0)
    eps_qv = -fracent*(qv - q0)
    eps_qt = -fracent*(qt - q0)-prate*(qt-qv)
    term1 = -B
    
    term2 = - Q_M*(Lv + Li*OMEGA)*cr['G']/(Rm0*T0)
    
    term3 = -cr['G']
    term4 = (cpm - Li*(qt - qv)*dOMEGA)*eps_T
    term5 = (Lv + Li*OMEGA)*(eps_qv + (qv/(1-qt))*eps_qt)

    term6 = cpm
    term7 = -Li*(qt - qv)*dOMEGA
    term8 = (Lv + Li*OMEGA)*(-dOMEGA*(qvv - qvi) + (1/(cr['Rv']*(T**2)))*(L_M))
    gamma_m = ( term1 + term2 + term3 + term4 + term5)/(term6 + term7 + term8)
    return gamma_m


# ----------------------------------------------------------------------------------------------------------------------------

def lift_parcel_adiabatic(T0, p0, q0, start_loc, fracent, prate, z0, T1, T2):
    """
    Compute lifted parcel properties using the unsaturated and saturated lapse
    rate formulas from (Peters et al. 2022).
    https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 
    
    Parameters
    ----------
    T0 : array_like
        Sounding profile of temperature (in K)
    p0 : array_like
        Sounding profile of pressure (in Pa)
    q0 : array_like
        Sounding profile of water vapor mass fraction (specific humidity), in kg/kg
    start_loc : int
        Index of the parcel starting location (set to 1 for the lowest level in
        the sounding)
    fracent : float
        Fractional entrainment rate (in m^-1)
    prate : float
        Precipitation rate (in m^-1) large values make parcel more pseudoadiabatic, 
        small values make parcel more adiabatic.  I usually just set it to 0 to
        get an adiabatic parcel
    z0 : array_like
        Sounding profile of height above ground level (AGL, first height should be 0 m)
    T1 : float
        Warmest mixed-phase temperature
    T2 : float
        Coldest mixed-phase temperature
    
    Returns
    -------
    T_lif : array_like
        Lifted parcel temperature
    Qv_lif : array_like
        Lifted parcel water vapor mass fraction (specific humidity), in kg/kg
    Qt_lif : array_like
        Lifted parcel total water mass fraction (cloud water), in kg/kg
    B_lif : array_like
        Lifted parcel buoyancy, computed using Eq. B6 in (Peters et al. 2022)
        (accounts for virtual temperature and loading effects)
    """
        
    # compute the moist static energy (MSE)
    MSE = compute_moist_static_energy(T0, q0, z0)
    # find the the index of the height of minimum MSE
    mn_hgt = np.where(MSE==np.min(MSE))
    # mn_hgt = int(mn_hgt[0])
    
    T_lif = np.zeros(T0.shape)*np.nan   # temperature of the lifted parcel
    Qv_lif = np.zeros(T0.shape)*np.nan  # water vapor mass fraction of the lifted parcel (specific humidity), in kg/kg
    Qt_lif = np.ones(T0.shape)*np.nan   # total water mass fraction of the lifted parcel

    # set initial values to that of the environment
    if start_loc>0:
        T_lif[0:start_loc+1] = T0[0:start_loc+1]
        Qv_lif[0:start_loc+1] = q0[0:start_loc+1]
        Qt_lif[0:start_loc+1] = Qv_lif[0:start_loc+1]
    else:
        T_lif[0] = T0[0]
        Qv_lif[0] = q0[0]
        Qt_lif[0] = Qv_lif[0]


    q_sat_prev = 0
    B_run = 0
    iz = start_loc

    #
    #for iz in np.arange(start_loc+1, z0.shape[0]):
    #
    #
    # I REVISED THIS A BIT.  TO MAKE THE CODE FASTER, I HAVE THE CALCULATION CUT OUT WHEN THE INTEGRATED NEGATIVE BUOYANCY ("BRUN") 
    # BECOMES MORE NEGATIVE THAN THE TOTAL INTEGRATED POSITIVE BUOYANCY.  I RESTRICT THIS TO ONLY HAPPEN AFTER WE HAVE PASSED 
    # THE HEIGHT OF MINIMUM MSE.  UNCOMMENT THE FOR LOOP ABOVE AND COMMENT OUT THE WHILE LOOP IF YOU JUST WANT TO INTEGRATE TO THE TOP OF THE SOUNDING.
    # THE +25 PART IN THE WHILE STATEMENT IS A PAD ON B_RUN (THE NEGATIVE CAPE HAS TO BE 25 J/KG LESS THAN THE POSITIVE CAPE TO KILL THE LOOP)
    # while iz<(z0.shape[0])-1 and (z0[iz]<z0[mn_hgt] or (B_run+25)>0):
    # while iz<(z0.shape[0])-1:
    while iz < (z0.shape[0])-1 and (z0[iz] < z0[mn_hgt] or (B_run+250) > 0):
        iz = iz + 1
        q_sat=(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1], p0[iz-1], T1, T2, 1)
        if Qv_lif[iz-1]<q_sat: # if we are unsaturated, go up at the unsaturated adiabatic lapse rate (eq. 19 in Peters et al. 2022)
            
        
        
            T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1]) * drylift(T_lif[iz-1], Qv_lif[iz-1], T0[iz-1], q0[iz-1], fracent)
            Qv_lif[iz] = Qv_lif[iz-1] - (z0[iz] - z0[iz-1]) * fracent * (Qv_lif[iz-1] - q0[iz-1])
            Qt_lif[iz] = Qv_lif[iz]
            q_sat = (1 - Qt_lif[iz]) * compute_rsat(T_lif[iz], p0[iz], T1, T2, 1)
            
            if Qv_lif[iz] >= q_sat: # if we hit saturation, split the vertical step into two stages.  The first stage advances at the saturated lapse rate to the saturation point, and the second stage completes the grid step at the moist lapse rate

                satrat = (Qv_lif[iz] - q_sat_prev)/(q_sat - q_sat_prev)
                dz_dry = satrat * (z0[iz] - z0[iz-1])
                dz_wet = (1-satrat) * (z0[iz] - z0[iz-1])


                
                T_halfstep = T_lif[iz-1] + dz_dry * drylift(T_lif[iz-1], Qv_lif[iz-1], T0[iz-1], q0[iz-1], fracent)
                Qv_halfstep = Qv_lif[iz-1] - dz_dry * fracent * (Qv_lif[iz-1] - q0[iz-1])
                Qt_halfstep = Qv_lif[iz]
                p_halfstep = p0[iz-1]*satrat + p0[iz] * (1-satrat)
                T0_halfstep = T0[iz-1]*satrat + T0[iz] * (1-satrat)
                Q0_halfstep = q0[iz-1]*satrat + q0[iz] * (1-satrat)

                T_lif[iz] = T_halfstep + dz_wet \
                            * moislif(T_halfstep, Qv_halfstep, (1-Qt_halfstep) \
                            *compute_rsat(T_halfstep, p_halfstep, T1, T2, 0), (1-Qt_halfstep) \
                            *compute_rsat(T_halfstep, p_halfstep, T1, T2, 2), p_halfstep, T0_halfstep, Q0_halfstep, Qt_halfstep, fracent, prate, T1, T2)
                
                
                Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1]) * fracent * (Qt_halfstep - Q0_halfstep)
                Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz], p0[iz], T1, T2, 1)

                if Qt_lif[iz] < Qv_lif[iz]:
                    Qv_lif[iz] = Qt_lif[iz]

            q_sat_prev = q_sat
            
        else: # if we are already at saturation, just advance upwacr['Rd'] using the saturated lapse rate (eq. 24 in Peters et al. 2022)
            T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1]) \
                        * moislif(T_lif[iz-1], Qv_lif[iz-1], (1-Qt_lif[iz-1]) \
                        * compute_rsat(T_lif[iz-1], p0[iz-1], T1, T2, 0), (1-Qt_lif[iz-1])\
                        * compute_rsat(T_lif[iz-1], p0[iz-1], T1, T2, 2), p0[iz-1], T0[iz-1], q0[iz-1], Qt_lif[iz-1], fracent, prate, T1, T2)
                     
             
            Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
            Qv_lif[iz] = (1-Qt_lif[iz]) * compute_rsat(T_lif[iz], p0[iz], T1, T2, 1)
            
            if Qt_lif[iz] < Qv_lif[iz]:
                Qv_lif[iz] = Qt_lif[iz]

        B_run = B_run + (cr['G']*T_lif[iz]*(1 + (cr['Rv']/cr['Rd'])*Qv_lif[iz] - Qt_lif[iz])/(T0[iz]*(1 + (cr['Rv']/cr['Rd'])*q0[iz] - q0[iz])) - cr['G'])*(z0[iz]-z0[iz-1])

    T_rho_lif = T_lif*(1 + (cr['Rv']/cr['Rd'])*Qv_lif - Qt_lif)
    T_0_lif = T0*(1 + (cr['Rv']/cr['Rd'] - 1)*q0)
    # T_rho_lif=T_lif*(1 - Qt_lif + Qv_lif)/( 1 + (epsilon - 1)/( ( epsilon*(1 - Qt_lif)/Qv_lif - 1) ) )
    # T_0_lif=T0/( 1 + (epsilon - 1)/( ( epsilon*(1 - q0)/q0 - 1) ) )
    
    B_lif = cr['G'] * (T_rho_lif - T_0_lif)/T_0_lif
    
    
    return T_lif, Qv_lif, Qt_lif, B_lif

# ----------------------------------------------------------------------------------------------------------------------------


def compute_CAPE_AND_CIN(T0, p0, q0, start_loc, fracent, prate, z0, T1, T2):

    """
    This function computes CAPE and CIN, along with the height of the Lifted
    Condensation Level (LCL) and the Equilibrium Level (EL).

    Parameters
    ----------
    T0 : float
        Temperature profile in K
    p0 : float
        Pressure profile in Pa
    q0 : float
        Water vapor mass fraction profile (specific humidity) in kg/kg
    start_loc : int
        Index of the parcel starting location
    fracent : float
        Fractional entrainment rate in m^-1
    prate : float
        Precipitation rate in m^-1
    z0 : float
        Height of the profile in m
    T1 : float
        Warmest mixed-phase temperature
    T2 : float
        Coldest mixed-phase temperature

    Returns
    -------
    CAPE : float
        Convective Available Potential Energy in J/kg
    CIN : float
        Convective Inhibition in J/kg
    LFC : float
        Lifted Condensation Level in m
    EL : float
        Equilibrium Level in m

    Notes
    -----
    If the parcel does not have any positive buoyancy, CAPE and CIN will be
    zero and the LCL and EL will be set to NaN.
    """

    # compute lifted parcel buoyancy
    _, _, _, B_lifted_parcel = lift_parcel_adiabatic(T0, p0, q0, start_loc, fracent, prate, z0, T1, T2)
    
    if np.nanmax(B_lifted_parcel)>0:

        # CAPE will be the total integrated positive buoyancy
        B_pos = np.zeros(B_lifted_parcel.shape)
        B_pos[:] = B_lifted_parcel[:]
        B_pos[np.where(B_pos<0)]=0
        dz = z0[1:z0.shape[0]] - z0[0:z0.shape[0]-1]
        CAPE = np.nansum(0.5*B_pos[0:z0.shape[0]-1]*dz + 0.5*B_pos[1:z0.shape[0]]*dz)
        
        # CIN will be the total negative buoyancy below the height of maximum
        B_neg = np.zeros(B_lifted_parcel.shape)
        B_neg[:] = B_lifted_parcel[:]
        mx = np.nanmax(B_lifted_parcel)
        imx = np.where(B_lifted_parcel==mx)
        imx=imx[0][0]
        B_neg[0:imx]=np.minimum( B_neg[0:imx], 0 )
        B_neg[imx:z0.shape[0]]= 0
        CIN = np.nansum(0.5*B_neg[0:z0.shape[0]-1]*dz + 0.5*B_neg[1:z0.shape[0]]*dz)
        
        # LFC will be the last instance of negative buoyancy before the
        # continuous integral that contains the maximum in buoyancy
        fneg = np.where(B_lifted_parcel<0)
        fneg=fneg[0]
        inn = np.where(fneg<imx)
        inn = inn[0]
        fneg = fneg[inn]
        if len(fneg)>0:
            LFC = 0.5*z0[np.max(fneg)] + 0.5*z0[np.max(fneg)+1]
        else:
            LFC = z0[start_loc]
        
        # EL will be last instance of positive buoyancy
        fpos = np.where(B_lifted_parcel>0)
        fpos=fpos[0]
        EL = 0.5*z0[np.max(fpos)] + 0.5*z0[np.max(fpos)+1]
    else:
        CAPE = 0
        CIN = 0
        LFC = np.nan
        EL = np.nan

    return CAPE, CIN, LFC, EL


# ----------------------------------------------------------------------------------------------------------------------------
# compute NCAPE
def compute_NCAPE(T0, p0, q0, z0, T1, T2, LFC, EL):
    """
    Compute NCAPE (Normalized CAPE) and other associated values.

    Parameters
    ----------
    T0 : numpy.ndarray
        Temperature of the parcel in K
    p0 : numpy.ndarray
        Pressure of the parcel in Pa
    q0 : numpy.ndarray
        Total water mass fraction of the parcel in kg/kg
    z0 : numpy.ndarray
        Height above ground level of the parcel in m
    T1 : float
        Warmest mixed-phase temperature
    T2 : float
        Coldest mixed-phase temperature
    LFC : float
        Height of the level of free convection in m
    EL : float
        Height of the equilibrium level in m

    Returns
    -------
    NCAPE : float
        Normalized CAPE in J/kg
    MSE0_star : numpy.ndarray
        Saturated moist static energy in J/kg
    MSE0bar : numpy.ndarray
        Average moist static energy in J/kg
    """
    if np.isnan(LFC) or np.isnan(EL):
        return np.nan, np.nan, np.nan

    # compute the moist static energy (MSE)
    MSE0 = compute_moist_static_energy(T0, q0, z0)
    
    # compute the saturated moist static energy
    rsat = compute_rsat(T0, p0, T1, T2, 0)
    qsat = (1 - rsat)*rsat
    MSE0_star = compute_moist_static_energy(T0, qsat, z0)
    
    # compute MSE0_BAR
    MSE0bar = np.zeros(MSE0.shape)
    #for iz in np.arange(0, MSE0bar.shape[0], 1):
     #   MSE0bar[iz]=np.mean(MSE0[1:iz])
        
    MSE0bar[0] = MSE0[0]
    for iz in np.arange(1, MSE0bar.shape[0], 1):
        MSE0bar[iz] = 0.5*np.sum((MSE0[0:iz] + MSE0[1:iz+1])*(z0[1:iz+1]-z0[0:iz]))/(z0[iz]-z0[0])
    
    int_arg = - (cr['G']/(cr['cpl']*T0)) * (MSE0bar - MSE0_star)
    ddiff = abs(z0-LFC)
    mn = np.min(ddiff)
    ind_LFC = np.where(ddiff == mn)[0][0]
    ddiff = abs(z0-EL)
    mn = np.min(ddiff)
    ind_EL = np.where(ddiff == mn)[0][0]

    
    NCAPE = np.maximum(np.nansum( (0.5*int_arg[ind_LFC:ind_EL-1] +
                                   0.5*int_arg[ind_LFC+1:ind_EL])*
                                  (z0[ind_LFC+1:ind_EL] - z0[ind_LFC:ind_EL-1]) ), 0)

    return NCAPE, MSE0_star, MSE0bar

# ----------------------------------------------------------------------------------------------------------------------------


def compute_VSR(z0, u0, v0):
    # compute 0-1 km storm-relative flow (V_SR) using the storm motion
    # estimate of Bunkers et al. (2000)
    # https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2
    
    f6000 = np.where(z0<=6000)[0]
    meanx=np.mean(u0[f6000])
    meany=np.mean(v0[f6000])
    
    f0500 = np.where(z0<=500)[0]
    lowx=np.mean(u0[f0500])
    lowy=np.mean(v0[f0500])
    
    f560 = np.where(np.logical_and(z0<=6000, z0>=5500))[0]
    highx=np.mean(u0[f560])
    highy=np.mean(v0[f560])
    BK_SHRx=highx-lowx
    BK_SHRy=highy-lowy
    BK_mag=np.sqrt(BK_SHRx**2 + BK_SHRy**2)
    BK_dirx=BK_SHRx/BK_mag
    BK_diry=BK_SHRy/BK_mag
    BK_orthx=BK_diry*7.5
    BK_orthy=-BK_dirx*7.5


    SR_mean_u = u0 - meanx
    SR_mean_v = v0 - meany
    dudz=np.zeros(u0.shape)
    dvdz=np.zeros(v0.shape)
    dudz[1:dudz.shape[0]-1] = (u0[2:dudz.shape[0]]-u0[0:dudz.shape[0]-2])/(z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2])
    dudz[0] = 2*dudz[1]-dudz[2]
    dvdz[1:dudz.shape[0]-1] = (v0[2:dudz.shape[0]]-v0[0:dudz.shape[0]-2])/(z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2])
    dvdz[0] = 2*dvdz[1]-dvdz[2]
    f1000 = np.where(z0<=1000)[0]
    SRH_mean = abs(np.mean(-SR_mean_u[f1000]*dvdz[f1000] + SR_mean_v[f1000]*dudz[f1000])*1000.0)
    
    
    propfac = min(SRH_mean/150, 1)
    propfac = 1


    C_x = meanx+propfac*BK_orthx
    C_y = meany+propfac*BK_orthy
    
    u_sr = u0 - C_x
    v_sr = v0 - C_y
    
    f1000 = np.where(z0<=1000)[0]
    V_SR = np.nanmean(np.sqrt(u_sr[f1000]**2 + v_sr[f1000]**2))
    return V_SR, C_x, C_y


def compute_VSR_DIFF(z0, u0, v0, rho0, EL, B_pos):
    # compute 0-1 km storm-relative flow (V_SR) using the storm motion
    # estimate of Bunkers et al. (2000)
    # https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2
    
    zdiff = (z0 - EL)**2
    ind_top = np.where(zdiff==np.min(zdiff))[0][0]
    inds_avg = np.arange(0, ind_top, 1)
    
    meanx = np.nanmean(B_pos[inds_avg]*rho0[inds_avg]*u0[inds_avg])/np.nanmean(B_pos[inds_avg]*rho0[inds_avg])
    meany = np.nanmean(B_pos[inds_avg]*rho0[inds_avg]*v0[inds_avg])/np.nanmean(B_pos[inds_avg]*rho0[inds_avg])
    
    #meanx = np.nanmean(rho0[inds_avg]*u0[inds_avg])/np.nanmean(rho0[inds_avg])
    #meany = np.nanmean(rho0[inds_avg]*v0[inds_avg])/np.nanmean(rho0[inds_avg])
    
    f6000 = np.where(z0<=6000)[0]
    #meanx=np.mean(u0[f6000])
    #meany=np.mean(v0[f6000])
    
    f0500 = np.where(z0<=500)[0]
    lowx=np.mean(u0[f0500])
    lowy=np.mean(v0[f0500])
    
    f560 = np.where(np.logical_and(z0<=6000, z0>=5500))[0]
    highx=np.mean(u0[f560])
    highy=np.mean(v0[f560])
    BK_SHRx=highx-lowx
    BK_SHRy=highy-lowy
    BK_mag=np.sqrt(BK_SHRx**2 + BK_SHRy**2)
    BK_dirx=BK_SHRx/BK_mag
    BK_diry=BK_SHRy/BK_mag
    BK_orthx=BK_diry*7.5
    BK_orthy=-BK_dirx*7.5


    SR_mean_u= u0 - meanx
    SR_mean_v= v0 - meany
    dudz=np.zeros(u0.shape)
    dvdz=np.zeros(v0.shape)
    dudz[1:dudz.shape[0]-1]= ( u0[2:dudz.shape[0]]-u0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dudz[0]=2*dudz[1]-dudz[2]
    dvdz[1:dudz.shape[0]-1]= ( v0[2:dudz.shape[0]]-v0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dvdz[0]=2*dvdz[1]-dvdz[2]
    f1000 = np.where(z0<=1000)[0]
    SRH_mean = abs(np.mean(-SR_mean_u[f1000]*dvdz[f1000] + SR_mean_v[f1000]*dudz[f1000])*1000.0)
    
    
    propfac=min(SRH_mean/150, 1)


    C_x=meanx+propfac*BK_orthx
    C_y=meany+propfac*BK_orthy
    
    u_sr = u0 - C_x
    v_sr = v0 - C_y
    
    f1000 = np.where(z0<=1000)[0]
    V_SR = np.nanmean(np.sqrt(  u_sr[f1000]**2 + v_sr[f1000]**2  ))
    return V_SR


# ----------------------------------------------------------------------------------------------------------------------------

def compute_ETILDE(CAPE, NCAPE, V_SR, EL, L):

    if np.isnan(CAPE) or np.isnan(NCAPE) or np.isnan(V_SR) or np.isnan(EL) or np.isnan(L):
        return np.nan, np.nan, np.nan
    # these are a bunch of constant paramters set for the ecape calculation
    H=EL
    sigma = 1.1
    alpha=0.8
    pitchfork = cr['ksq']*(alpha**2)*(np.pi**2)*L/(4*cr['Pr']*(sigma**2)*H)
    vsr_tilde = V_SR/np.sqrt(2*CAPE)
    N_tilde = NCAPE/CAPE
    
    # equation solves for the nondimensional ecape (e_tilde_a in the paper)
    E_tilde = vsr_tilde**2 + ( -1 - pitchfork - (pitchfork/(vsr_tilde**2 ))*N_tilde + \
                              np.sqrt((1 + pitchfork + (pitchfork/(vsr_tilde**2 ))*N_tilde)**2 + \
                                      (4*(pitchfork/(vsr_tilde**2 ))*(1 - pitchfork*N_tilde) ) ) )/( 2*pitchfork/(vsr_tilde**2) )
        
    E_tilde_ = E_tilde - vsr_tilde**2
        
    varepsilon = 2*((1 - E_tilde_)/(E_tilde_ + N_tilde))/(EL)  
    

    #eps = 2*cr['ksq']*L/(EL*cr['Pr'])
    
    #Rm2 = ( (alpha*np.pi/(sigma) )**2 )*( E_tilde/vsr_tilde + 1)
    #Radius =  EL*Rm2**(-1/2)
    #varepsilon = 2*cr['ksq']*L/(cr['Pr']*Radius**2 )
    
    #Radius=Radius/2
    
    # Fractional Entrainment Rate
    # varepsilon = 0.65*eps*(alpha**2)*(np.pi**2)*E_tilde/(4*(sigma**2)*EL*(vsr_tilde**2 ) )

    # Radius
    Radius = np.sqrt(2*cr['ksq']*L/(cr['Pr']*varepsilon))

    return E_tilde, varepsilon, Radius


