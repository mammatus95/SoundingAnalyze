# ----------------------------------------------------------------------------------------------------------------------------
#                 Run ECAPE
# Author : Morten Kretschmer
# 
# ----------------------------------------------------------------------------------------------------------------------------

import numpy as np
import time
# ----------------------------------------------------------------------------------------------------------------------------

from src.meteolib import cr, temp_at_mixrat, q_to_mixrat
from src.cm1_lib import read_modelsounding, interpolating_sounding, calculate_density_potential_temperature
from src.cm1_lib import calculate_PII, calculate_pressure, calculate_temperature_density, calculate_density, calculate_temperature
from src.ecape_lib import CI_model, compute_CAPE_AND_CIN, compute_w, compute_NCAPE, compute_VSR, compute_ETILDE
from src.ecape_lib import lift_parcel_adiabatic
from src.plotlib import plot_stuve_cm1

# ----------------------------------------------------------------------------------------------------------------------------


def main():


    sound_filename = 'src/example/sounding_base.txt'
    #sound_filename = 'src/example/sounding_noshear.txt'

    z_env, Th_env, qv_env, u_env, v_env, p_sfc = read_modelsounding(sound_filename)
    z_env, Th_env, qv_env, u_env, v_env = interpolating_sounding(z_env, Th_env, qv_env, u_env, v_env)
    mr_env = q_to_mixrat(qv_env)

    Th_rho = calculate_density_potential_temperature(Th_env, qv_env)

    PII = calculate_PII(z_env, Th_rho, p_sfc*100.0)

    pres_env = calculate_pressure(PII)

    # T_rho = calculate_temperature_density(PII, Th_env, qv_env)

    # rho = calculate_density(pres_env, T_rho)

    T_env = calculate_temperature(pres_env, Th_env)

    TD_env = temp_at_mixrat(mr_env*1000, pres_env/100) # degree C

    T1 = 273.15
    T2 = 253.15

    # time profiling ECAPE
    start = time.time()

    # GET THE SURFACE-BASED CAPE, CIN, LFC, EL
    CAPE, CIN, LFC, EL = compute_CAPE_AND_CIN(T_env, pres_env, qv_env, 0, 0, 0, z_env, T1, T2)
    
    # GET NCAPE, WHICH IS NEEDED FOR ECAPE CALULATION
    NCAPE, MSE0_star, MSE0bar = compute_NCAPE(T_env, pres_env, qv_env, z_env, T1, T2, LFC, EL)

    # GET THE 0-1 KM MEAN STORM-RELATIVE WIND, ESTIMATED USING BUNKERS METHOD FOR RIGHT-MOVER STORM MOTION
    V_SR, C_x, C_y = compute_VSR(z_env, u_env, v_env)

    # GET E_TILDE, WHICH IS THE RATIO OF ECAPE TO CAPE.  ALSO, VAREPSILON IS THE FRACITONAL ENTRAINMENT RATE, AND RADIUS IS THE THEORETICAL UPRAFT RADIUS
    E_tilde, varepsilon, Radius = compute_ETILDE(CAPE, NCAPE, V_SR, EL, 250)

    fracent=varepsilon
    #prate=3e-5
    ECAPE, ECIN, ELFC, EEL = compute_CAPE_AND_CIN(T_env, pres_env, qv_env, 0, fracent, 0, z_env, T1, T2)
    
    # end time profiling ECAPE
    end = time.time()

    # print results
    print(f"CAPE: {CAPE} CIN: {CIN} LFC: {LFC} EL: {EL}")
    print(f"ECAPE: {ECAPE} ECIN: {ECIN} ELFC: {ELFC} EEL: {EEL}")
    print(f"time elapsed: {(end - start):.2f} sec") # in seconds


    WCAPE, WCIN, WLFC, WEL = compute_w(T_env, pres_env, qv_env, 0, 0, 0, z_env, T1, T2, 1000, u_env, v_env, V_SR)

    # CI THEORETICAL MODEL
    R_TS = CI_model(T_env, pres_env, qv_env, z_env, u_env, v_env, T1, T2, np.arange(300, 6000, 100), 20, 250, 0)





    # GET THE ECAPE AND ECAPE WITH ADIABATIC PARCEL LIFT
    print(f"CAPE: {CAPE} WCAPE: {WCAPE} NCAPE: {NCAPE} E_tilde: {E_tilde} ECAPE: {ECAPE}")

    # calculate parcel trajectories for plotting
    T_lif, Qv_lif, Qt_lif, _ = lift_parcel_adiabatic(T_env, pres_env, qv_env, 0, 0, 0, z_env, T1, T2)
    T_lif_ecape, Qv_lif_ecape, Qt_lif_ecape, _ = lift_parcel_adiabatic(T_env, pres_env, qv_env, 0, fracent, 0, z_env, T1, T2)

    plot_stuve_cm1(T_env - cr['ZEROCNK'], TD_env, pres_env, T_lif, Qv_lif, Qt_lif, T_lif_ecape, Qv_lif_ecape, Qt_lif_ecape)

if __name__ == '__main__':
    print("Running ECAPE")
    main()