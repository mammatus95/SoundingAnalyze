# ----------------------------------------------------------------------------------------------------------------------------
#                 Run ECAPE
# Author : Morten Kretschmer
# 
# ----------------------------------------------------------------------------------------------------------------------------

import time
# ----------------------------------------------------------------------------------------------------------------------------

from src.meteolib import cr, temp_at_mixrat, q_to_mixrat
from src.cm1_lib import read_modelsounding, interpolating_sounding, calculate_density_potential_temperature
from src.cm1_lib import calculate_PII, calculate_pressure, calculate_temperature  # calculate_temperature_density, calculate_density
from src.ecape_lib import compute_CAPE_AND_CIN, compute_NCAPE, compute_VSR, compute_ETILDE
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

    # get surface based CAPE, CIN, LFC, EL
    CAPE, CIN, LFC, EL = compute_CAPE_AND_CIN(T_env, pres_env, qv_env, 0, 0, 0, z_env, T1, T2)
    
    # get NCAPE, which is needed for ECAPE calulation.
    NCAPE, MSE0_star, MSE0bar = compute_NCAPE(T_env, pres_env, qv_env, z_env, T1, T2, LFC, EL)

    # get the 0-1 km mean storm-relative wind, estimated using bunkers method for right-mover storm motion
    V_SR, _, _ = compute_VSR(z_env, u_env, v_env)

    # get E_TILDE, which is the ratio of ecape to cape.
    # Also, varepsilon is the fracitonal entrainment rate, and radius is the theoretical upraft radius
    E_tilde, varepsilon, Radius = compute_ETILDE(CAPE, NCAPE, V_SR, EL, 250)

    fracent=varepsilon
    # prate=3e-5
    ECAPE, ECIN, ELFC, EEL = compute_CAPE_AND_CIN(T_env, pres_env, qv_env, 0, fracent, 0, z_env, T1, T2)
    
    # end time profiling ECAPE
    end = time.time()

    # print results
    print(f"CAPE: {CAPE} CIN: {CIN} LFC: {LFC} EL: {EL}")
    print(f"ECAPE: {ECAPE} ECIN: {ECIN} ELFC: {ELFC} EEL: {EEL}")
    print(f"time elapsed: {(end - start)*1000:.3f} ms") # in milliseconds
    #exit(0)


    # calculate parcel trajectories for plotting
    T_lif, Qv_lif, Qt_lif, _ = lift_parcel_adiabatic(T_env, pres_env, qv_env, 0, 0, 0, z_env, T1, T2)
    T_lif_ecape, Qv_lif_ecape, Qt_lif_ecape, _ = lift_parcel_adiabatic(T_env, pres_env, qv_env, 0, fracent, 0, z_env, T1, T2)

    # plotting
    plot_stuve_cm1(T_env - cr['ZEROCNK'], TD_env, pres_env, T_lif, Qv_lif, Qt_lif, T_lif_ecape, Qv_lif_ecape, Qt_lif_ecape)

if __name__ == '__main__':
    print("Running ECAPE")
    main()