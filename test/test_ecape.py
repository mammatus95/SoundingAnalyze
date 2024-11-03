#!/usr/bin/python3
import unittest
import numpy as np
# project modul
from src.ecape_lib import omega, domega, get_qs, compute_moist_static_energy, compute_rsat, moislif
from src.ecape_lib import compute_LCL, compute_LCL_NUMERICAL
from src.meteolib import cr


# ----------------------------------------------------------------------------------------------------------------------------


class TestOmega(unittest.TestCase):
    def setUp(self):
        self.T1 = 273.15
        self.T2 = 253.15

    def test_omega_1(self):
        T = 263.15
        expected_result = 0.5
        self.assertAlmostEqual(omega(T, self.T1, self.T2), expected_result)

    def test_omega_2(self):
        T = 268.15
        expected_result = 0.25
        self.assertAlmostEqual(omega(T, self.T1, self.T2), expected_result)

    def test_omega_3(self):
        T = 258.15
        expected_result = 0.75
        self.assertAlmostEqual(omega(T, self.T1, self.T2), expected_result)

    def test_omega_4(self):
        T = 253.15
        expected_result = 1.0
        self.assertAlmostEqual(omega(T, self.T1, self.T2), expected_result)

    def test_omega_5(self):
        T = 303.15
        expected_result = 0.0
        self.assertAlmostEqual(omega(T, self.T1, self.T2), expected_result)

# ----------------------------------------------------------------------------------------------------------------------------


class TestDomegaFunction(unittest.TestCase):
    def setUp(self):
        self.T1 = 273.15
        self.T2 = 253.15

    def test_domega_within_range(self):
        T = 270.0
        expected_result = -0.05
        self.assertAlmostEqual(domega(T, self.T1, self.T2), expected_result)

    def test_domega_less_than_T2(self):
        T = 240.0
        expected_result = 0.0
        self.assertAlmostEqual(domega(T, self.T1, self.T2), expected_result)

    def test_domega_greater_than_T1(self):
        T = 300.0
        expected_result = 0.0
        self.assertAlmostEqual(domega(T, self.T1, self.T2), expected_result)

    def test_domega_T1_equal_T2(self):
        T = 260.5
        T1 = 273.15
        T2 = 273.15
        with self.assertRaises(ZeroDivisionError):
            domega(T, T1, T2)
        
    def test_realistic_domega(self):
        T = 260.0
        expected_result = -0.05
        self.assertAlmostEqual(domega(T, self.T1, self.T2), expected_result)

# ----------------------------------------------------------------------------------------------------------------------------


class TestGetQsFunction(unittest.TestCase):

    def test_single_value_inputs(self):
        qt = 0.5
        rs = 0.2
        expected_result = 0.1
        self.assertAlmostEqual(get_qs(qt, rs), expected_result)

    def test_array_inputs(self):
        qt = np.array([0.5, 0.7])
        rs = np.array([0.2, 0.3])
        expected_result = np.array([0.1, 0.09])
        self.assertTrue(type(get_qs(qt, rs)), type(np.array([1,2])))
        np.testing.assert_array_almost_equal(get_qs(qt, rs), expected_result)

    def test_zero_input_qt(self):
        qt = 0
        rs = 0.2
        expected_result = 0.2
        self.assertAlmostEqual(get_qs(qt, rs), expected_result)

    def test_zero_input_rs(self):
        qt = 0.5
        rs = 0
        expected_result = 0
        self.assertAlmostEqual(get_qs(qt, rs), expected_result)

    def test_negative_input_qt(self):
        qt = -0.5
        rs = 0.2
        expected_result = 0.3
        self.assertAlmostEqual(get_qs(qt, rs), expected_result)

    def test_negative_input_rs(self):
        qt = 0.5
        rs = -0.2
        expected_result = -0.1
        self.assertAlmostEqual(get_qs(qt, rs), expected_result)

# ----------------------------------------------------------------------------------------------------------------------------


class TestComputeMoistStaticEnergy(unittest.TestCase):
    def setUp(self):
        global cr
        cr = {'cpl': 4218.0, 'xlv': 2501000, 'G': 9.80665}  # example values for cr
        self.cr = cr


    def test_MSE(self):
        T_env = 300.0
        mr_env = 0.5
        z0 = 100.0
        expected_result = cr['cpl']*T_env + cr['xlv']*mr_env + cr['G']*z0
        self.assertAlmostEqual(compute_moist_static_energy(T_env, mr_env, z0), expected_result)

# ----------------------------------------------------------------------------------------------------------------------------


class TestComputeRsat(unittest.TestCase):
    def setUp(self):
        global cr
        # cpl 4190.0 4218.0 4186.0 J/kg/K   liquid water specific heat capacity
        cr  = {'Rd': 287.05, 'Rv': 461.51, 'cpv': 1875.0, 'cpl': 4186.0, 'cpi': 1840.0, 'xlv': 2.5e6, 'xls': 2.834e6, 'ttrip': 273.16, 'eref': 611.73, 'G': 9.81}
        self.cr = cr
        self.T1 = 273.15
        self.T2 = 253.15

    def test_compute_rsat_liquid(self):
        T = 300.0
        p = 100000.0
        iceflag = 0
        expected_result = 0.022780664040943
        qsat = compute_rsat(T, p, self.T1, self.T2, iceflag)
        self.assertGreater(qsat, 0.0)
        self.assertAlmostEqual(qsat, expected_result)

    def test_compute_rsat_linear_combination(self):
        T = 300.0
        p = 100000.0
        iceflag = 1
        qsat = compute_rsat(T, p, self.T1, self.T2, iceflag)
        self.assertGreater(qsat, 0.0)

    def test_compute_rsat_ice_only(self):
        T = 250.0
        p = 100000.0
        iceflag = 2
        qsat = compute_rsat(T, p, self.T1, self.T2, iceflag)
        self.assertGreater(qsat, 0.0)

    def test_compute_rsat_invalid_iceflag(self):
        T = 300.0
        p = 100000.0
        iceflag = 3
        with self.assertRaises(ValueError):
            compute_rsat(T, p, self.T1, self.T2, iceflag)

# ----------------------------------------------------------------------------------------------------------------------------


class TestMoislif(unittest.TestCase):
    def test_typical_input(self):
        T = 288.15  # K
        qv = 0.014  # kg/kg
        qvv = 0.014  # kg/kg
        qvi = 0.0  # kg/kg
        pres_env = 101325  # Pa
        T_env = 288.15  # K
        mr_env = 0.01  # kg/kg
        qt = 0.01  # kg/kg
        fracent = 0.001  # m^-1
        prate = 0.0  # m^-1
        T1 = 273.15  # K
        T2 = 233.15  # K

        gamma_m = moislif(T, qv, qvv, qvi, pres_env, T_env, mr_env, qt, fracent, prate, T1, T2)
        self.assertLess(gamma_m, 0.0)  # lapse rate should be positive

    def test_upperair_input(self):
        T = 260.0  # K
        qv = 0.001  # kg/kg
        qvv = 0.001  # kg/kg
        qvi = 0.001  # kg/kg
        pres_env = 70000  # Pa
        T_env = 256.15  # K
        mr_env = 0.001  # kg/kg
        qt = 0.01  # kg/kg
        fracent = 0.001  # m^-1
        prate = 0.0  # m^-1
        T1 = 273.15  # K
        T2 = 233.15  # K

        gamma_m = moislif(T, qv, qvv, qvi, pres_env, T_env, mr_env, qt, fracent, prate, T1, T2)
        self.assertLess(gamma_m, 0.0)  # lapse rate should be positive

    """
    def test_invalid_input(self):
        T = -100.0  # K
        qv = 0.01  # kg/kg
        qvv = 0.01  # kg/kg
        qvi = 0.0  # kg/kg
        pres_env = 101325  # Pa
        T_env = 288.15  # K
        mr_env = 0.01  # kg/kg
        qt = 0.01  # kg/kg
        fracent = 0.001  # m^-1
        prate = 0.0  # m^-1
        T1 = 273.15  # K
        T2 = 233.15  # K

        with self.assertRaises(ValueError):
            moislif(T, qv, qvv, qvi, pres_env, T_env, mr_env, qt, fracent, prate, T1, T2)
    """


# ----------------------------------------------------------------------------------------------------------------------------


class TestComputeLCL(unittest.TestCase):
    def test_typical_values(self):
        T = 288.15  # K
        qv = 0.01  # kg/kg
        p = 101325.0  # Pa
        Z_LCL, T_LCL, P_LCL = compute_LCL(T, qv, p)
        self.assertAlmostEqual(Z_LCL, 103.35, places=2)
        self.assertAlmostEqual(T_LCL, 287.15, places=2)
        self.assertAlmostEqual(P_LCL, 100096.26, places=2)

    def test_dry_condition(self):
        T = 288.15  # K
        qv = 0.002  # kg/kg
        p = 101325  # Pa
        Z_LCL, T_LCL, P_LCL = compute_LCL(T, qv, p)
        self.assertAlmostEqual(Z_LCL, 2871.59, places=2)
        self.assertAlmostEqual(T_LCL, 260.20, places=2)
        self.assertAlmostEqual(P_LCL, 70854.48, places=2)



class TestComputeLCLNumerical(unittest.TestCase):

    def test_typical_values(self):
        T = 288.15  # K
        qv = 0.01  # kg/kg
        p = 101325  # Pa
        dz = 10  # m
        Z_LCL = compute_LCL_NUMERICAL(T, qv, p, dz)
        self.assertGreater(Z_LCL, 0)
        self.assertAlmostEqual(Z_LCL, 110.0, places=2)

    def test_dry_condition(self):
        T = 288.15  # K
        qv = 0.002  # kg/kg
        p = 101325  # Pa
        dz = 10  # m
        Z_LCL = compute_LCL_NUMERICAL(T, qv, p, dz)
        self.assertGreater(Z_LCL, 0)
        self.assertAlmostEqual(Z_LCL, 2870.0, places=2)

    """
    def test_zero_vertical_resolution(self):
        T = 288.15  # K
        qv = 0.01  # kg/kg
        p = 101325  # Pa
        dz = 0  # m
        with self.assertRaises(ZeroDivisionError):
            compute_LCL_NUMERICAL(T, qv, p, dz)

    def test_negative_vertical_resolution(self):
        T = 288.15  # K
        qv = 0.01  # kg/kg
        p = 101325  # Pa
        dz = -10  # m
        with self.assertRaises(ValueError):
            compute_LCL_NUMERICAL(T, qv, p, dz)
    """


# ----------------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    unittest.main()
