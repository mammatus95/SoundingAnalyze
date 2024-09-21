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

        self.expected_cape = 3339.7
        self.expected_cin = -50.4
        self.expected_lfc = 1750.0
        self.expected_el = 11750.0

        self.expected_ncape = 0

        self.expected_ecape = 2250.1
        self.expected_ecin = -51.5
        self.expected_elfc = 1750.0
        self.expected_eel = 11150.0

        self.T1 = 273.15
        self.T2 = 253.15

        # read sound and prepare function call of different CAPE kinds
        try:
            Z, Th, r_v, u, v, p_sfc = read_modelsounding(self.sound_filename)
        except FileNotFoundError:
            Z, Th, r_v, u, v, p_sfc = read_modelsounding(self.sound_filename2)

        Z, Th, r_v,self.u, self.v = interpolating_sounding(Z, Th, r_v, u, v)
        Th_rho = calculate_density_potential_temperature(Th, r_v)

        PII = calculate_PII(Z, Th_rho, p_sfc*100.0)

        pres = calculate_pressure(PII)
        self.T0 = calculate_temperature(pres, Th)
        self.T_rho = calculate_temperature_density(PII, Th, r_v)

        self.q0 = (1 - r_v)*r_v
        self.p0 = pres
        self.z0 = Z


    def test_CAPE(self):
        CAPE, CIN, LFC, EL = compute_CAPE_AND_CIN(self.T0, self.p0, self.q0, 0, 0, 0, self.z0, self.T1, self.T2)

        self.assertAlmostEqual(CAPE, self.expected_cape, delta=self.delta)
        self.assertAlmostEqual(CIN, self.expected_cin, delta=self.delta)
        self.assertAlmostEqual(LFC, self.expected_lfc)
        self.assertAlmostEqual(EL, self.expected_el)

        self.LFC = LFC
        self.EL = EL


    def test_NCAPE(self):
        NCAPE, _, _ = compute_NCAPE(self.T0, self.p0, self.q0, self.z0, self.T1, self.T2, self.LFC, self.EL)

        self.assertAlmostEqual(NCAPE, self.expected_ncape, delta=self.delta)

        self.NCAPE = NCAPE


    def test_ECAPE(self):
        V_SR, _, _ = compute_VSR(self.z0, self.u, self.v)
        E_tilde, varepsilon, _ = compute_ETILDE(self.CAPE, self.NCAPE, V_SR, self.EL, 250)

        fracent=varepsilon
        ECAPE, ECIN, ELFC, EEL = compute_CAPE_AND_CIN(self.T0, self.p0, self.q0, 0, fracent, 0, self.z0, self.T1, self.T2)
        # ECAPE = E_tilde*self.CAPE
        self.assertAlmostEqual(E_tilde*self.CAPE, self.expected_ecape, delta=self.delta)
        self.assertAlmostEqual(ECAPE, self.expected_ecape, delta=self.delta)
        self.assertAlmostEqual(ECIN, self.expected_ecin, delta=self.delta)
        self.assertAlmostEqual(ELFC, self.expected_elfc)
        self.assertAlmostEqual(EEL, self.expected_eel)

# ----------------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    unittest.main()
