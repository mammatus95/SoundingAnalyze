# ----------------------------------------------------------------------------------------------------------------------------
#                 Run ECAPE
# Author : Morten Kretschmer
# 
# ----------------------------------------------------------------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------------------------------------------------------------

from src.meteolib import cr
from src.cm1_lib import read_modelsounding, interpolating_sounding, calculate_density_potential_temperature
from src.cm1_lib import calculate_PII, calculate_pressure, calculate_temperature_density, calculate_density, calculate_temperature
from src.ecape_lib import CI_model, compute_CAPE_AND_CIN, compute_w, compute_NCAPE, compute_VSR, compute_ETILDE
from src.ecape_lib import lift_parcel_adiabatic
# ----------------------------------------------------------------------------------------------------------------------------


def main():

    sound_filename = 'src/example/sounding_base.txt'
    Z, Th, r_v, u, v, p_sfc = read_modelsounding(sound_filename)
    Z, Th, r_v, u, v = interpolating_sounding(Z, Th, r_v, u, v)


    Th_rho = calculate_density_potential_temperature(Th, r_v)

    PII = calculate_PII(Z, Th_rho, p_sfc)

    pres = calculate_pressure(PII)

    T_rho = calculate_temperature_density(PII, Th, r_v)

    rho = calculate_density(pres, T_rho)

    T0 = calculate_temperature(pres, Th)

    T1 = 273.15
    T2 = 253.15

    q0 = (1 - r_v)*r_v
    p0 = pres
    z0 = Z

    # GET THE SURFACE-BASED CAPE, CIN, LFC, EL
    CAPE, CIN, LFC, EL = compute_CAPE_AND_CIN(T0, p0, q0, 0, 0, 0, z0, T1, T2)
    
    # GET NCAPE, WHICH IS NEEDED FOR ECAPE CALULATION
    NCAPE, MSE0_star, MSE0bar = compute_NCAPE(T0, p0, q0, z0, T1, T2, LFC, EL)

    # GET THE 0-1 KM MEAN STORM-RELATIVE WIND, ESTIMATED USING BUNKERS METHOD FOR RIGHT-MOVER STORM MOTION
    V_SR, C_x, C_y = compute_VSR(z0, u, v)

    WCAPE, WCIN, WLFC, WEL = compute_w(T0, p0, q0, 0, 0, 0, z0, T1, T2, 1000, u, v, V_SR)


    # GET E_TILDE, WHICH IS THE RATIO OF ECAPE TO CAPE.  ALSO, VAREPSILON IS THE FRACITONAL ENTRAINMENT RATE, AND RADIUS IS THE THEORETICAL UPRAFT RADIUS
    E_tilde, varepsilon, Radius = compute_ETILDE(CAPE, NCAPE, V_SR, EL, 250)

    # CI THEORETICAL MODEL
    R_TS = CI_model(T0, p0, q0, z0, u, v, T1, T2, np.arange(300, 6000, 100), 20, 250, 0)

    
    fracent=varepsilon
    #prate=3e-5
    T_lif, Qv_lif, Qt_lif, B_lif = lift_parcel_adiabatic(T0, p0, q0, 0, fracent, 0, z0, T1, T2)
    ECAPE, ECIN, ELFC, EEL = compute_CAPE_AND_CIN(T0, p0, q0, 0, fracent, 0, z0, T1, T2)

    print(f"WCAPE: {WCAPE} NCAPE: {NCAPE} E_tilde: {E_tilde} ECAPE: {ECAPE} CAPE: {CAPE}")

if __name__ == '__main__':
    print("Running ECAPE")
    main()