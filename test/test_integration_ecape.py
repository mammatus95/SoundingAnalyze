#!/usr/bin/python3
import unittest
import numpy as np
# project modul
from src.cm1_lib import read_modelsounding, interpolating_sounding, calculate_density_potential_temperature
from src.cm1_lib import calculate_PII, calculate_pressure, calculate_temperature_density, calculate_temperature
from src.ecape_lib import compute_CAPE_AND_CIN, compute_NCAPE, compute_VSR, compute_ETILDE


# ----------------------------------------------------------------------------------------------------------------------------


class TestCAPE(unittest.TestCase):
    def setUp(self):
        self.sound_filename = './sounding_base.txt'
        self.sound_filename2 = 'src/example/sounding_base.txt'

        # absolute allowed difference esimate and expected value
        self.delta = 10.0 # J/kg 

        self.expected_cape = 3531.94
        self.expected_cin = -50.4
        self.expected_lfc = 1650.0
        self.expected_el = 11850.0

        self.expected_ncape = 0

        self.expected_ecape = 2328.29
        self.expected_ecin = -51.5
        self.expected_elfc = 1750.0
        self.expected_eel = 11150.0

        self.T1 = 273.15
        self.T2 = 253.15

        # read sound and prepare function call of different CAPE kinds
        try:
            z_env, Th_env, qv_env, u_env, v_env, p_sfc = read_modelsounding(self.sound_filename)
        except FileNotFoundError:
            z_env, Th_env, qv_env, u_env, v_env, p_sfc = read_modelsounding(self.sound_filename2)

        z_env, Th_env, qv_env, self.u_env, self.v_env = interpolating_sounding(z_env, Th_env, qv_env, u_env, v_env)
        Th_rho = calculate_density_potential_temperature(Th_env, qv_env)

        PII = calculate_PII(z_env, Th_rho, p_sfc*100.0)

        self.p_env = calculate_pressure(PII)
        self.T_env = calculate_temperature(self.p_env, Th_env)
        self.T_rho = calculate_temperature_density(PII, Th_env, qv_env)
        self.qv_env = qv_env
        self.z_env = z_env

        self.CAPE = 3548.941689903609
        self.CIN = -50.42575575901518
        self.LFC = 1750.0
        self.EL = 11750.0
        self.NCAPE = 0.0
        self.ECAPE = np.nan
        self.ECIN = np.nan
        self.ELFC = np.nan
        self.EEL = np.nan


    def test_CAPE(self):
        CAPE, CIN, LFC, EL = compute_CAPE_AND_CIN(self.T_env, self.p_env, self.qv_env, 0, 0, 0, self.z_env, self.T1, self.T2)

        self.assertAlmostEqual(CAPE, self.expected_cape, delta=self.delta)
        self.assertAlmostEqual(CIN, self.expected_cin, delta=self.delta)
        self.assertAlmostEqual(LFC, self.expected_lfc)
        self.assertAlmostEqual(EL, self.expected_el)

        self.CAPE = CAPE
        self.CIN = CIN
        self.LFC = LFC
        self.EL = EL


    def test_NCAPE(self):
        NCAPE, _, _ = compute_NCAPE(self.T_env, self.p_env, self.qv_env, self.z_env, self.T1, self.T2, self.LFC, self.EL)

        self.assertAlmostEqual(NCAPE, self.expected_ncape, delta=self.delta)

        self.NCAPE = NCAPE


    def test_ECAPE(self):
        V_SR, _, _ = compute_VSR(self.z_env, self.u_env, self.v_env)
        E_tilde, varepsilon, _ = compute_ETILDE(self.CAPE, self.NCAPE, V_SR, self.EL, 250)

        fracent=varepsilon
        ECAPE, ECIN, ELFC, EEL = compute_CAPE_AND_CIN(self.T_env, self.p_env, self.qv_env, 0, fracent, 0, self.z_env, self.T1, self.T2)
        # ECAPE = E_tilde*self.CAPE
        #self.assertAlmostEqual(E_tilde*self.CAPE, self.expected_ecape, delta=self.delta)
        self.assertAlmostEqual(ECAPE, self.expected_ecape, delta=self.delta)
        self.assertAlmostEqual(ECIN, self.expected_ecin, delta=self.delta)
        self.assertAlmostEqual(ELFC, self.expected_elfc)
        self.assertAlmostEqual(EEL, self.expected_eel)

# ----------------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    unittest.main()
