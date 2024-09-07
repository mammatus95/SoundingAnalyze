#!/usr/bin/python3
import unittest
import numpy as np
# project modul
from src.ecape_lib import omega, domega, get_qs, compute_moist_static_energy, compute_rsat, moislif
from src.meteolib import cr



class TestOmega(unittest.TestCase):
    def test_omega_1(self):
        T, T1, T2 = 0, -10, 10
        expected_result = 0.5
        self.assertAlmostEqual(omega(T, T1, T2), expected_result)

    def test_omega_2(self):
        T, T1, T2 = 5, -10, 10
        expected_result = 0.75
        self.assertAlmostEqual(omega(T, T1, T2), expected_result)

    def test_omega_3(self):
        T, T1, T2 = 15, -10, 10
        expected_result = 1
        self.assertAlmostEqual(omega(T, T1, T2), expected_result)

    def test_omega_4(self):
        T, T1, T2 = -5, -10, 10
        expected_result = 0.25
        self.assertAlmostEqual(omega(T, T1, T2), expected_result)


class TestDomegaFunction(unittest.TestCase):

    def test_domega_within_range(self):
        T = 0.5
        T1 = 0.0
        T2 = 1.0
        expected_result = -1.0
        self.assertAlmostEqual(domega(T, T1, T2), expected_result)

    def test_domega_less_than_T1(self):
        T = -1.0
        T1 = 0.0
        T2 = 1.0
        expected_result = 0.0
        self.assertAlmostEqual(domega(T, T1, T2), expected_result)

    def test_domega_greater_than_T2(self):
        T = 2.0
        T1 = 0.0
        T2 = 1.0
        expected_result = 0.0
        self.assertAlmostEqual(domega(T, T1, T2), expected_result)

    def test_domega_T1_equal_T2(self):
        T = 0.5
        T1 = 1.0
        T2 = 1.0
        with self.assertRaises(ZeroDivisionError):
            domega(T, T1, T2)


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



class TestComputeMoistStaticEnergy(unittest.TestCase):
    def setUp(self):
        global cr
        cr = {'cpl': 4218.0, 'xlv': 2501000, 'G': 9.80665}  # example values for cr
        self.cr = cr


    def test_MSE(self):
        T0 = 300.0
        q0 = 0.5
        z0 = 100.0
        expected_result = cr['cpl']*T0 + cr['xlv']*q0 + cr['G']*z0
        self.assertAlmostEqual(compute_moist_static_energy(T0, q0, z0), expected_result)



class TestComputeRsat(unittest.TestCase):
    def setUp(self):
        global cr
        # cpl 4190.0 4218.0 4186.0 J/kg/K   liquid water specific heat capacity
        cr  = {'Rd': 287.05, 'Rv': 461.51, 'cpv': 1875.0, 'cpl': 4186.0, 'cpi': 1840.0, 'xlv': 2.5e6, 'xls': 2.834e6, 'ttrip': 273.16, 'eref': 611.73, 'G': 9.81}
        self.cr = cr

    def test_compute_rsat_liquid(self):
        T = 300.0
        p = 100000.0
        T1 = 250.0
        T2 = 350.0
        iceflag = 0
        qsat = compute_rsat(T, p, T1, T2, iceflag)
        self.assertGreater(qsat, 0.0)

    def test_compute_rsat_linear_combination(self):
        T = 300.0
        p = 100000.0
        T1 = 250.0
        T2 = 350.0
        iceflag = 1
        qsat = compute_rsat(T, p, T1, T2, iceflag)
        self.assertGreater(qsat, 0.0)

    def test_compute_rsat_ice_only(self):
        T = 250.0
        p = 100000.0
        T1 = 200.0
        T2 = 300.0
        iceflag = 2
        qsat = compute_rsat(T, p, T1, T2, iceflag)
        self.assertGreater(qsat, 0.0)

    def test_compute_rsat_invalid_iceflag(self):
        T = 300.0
        p = 100000.0
        T1 = 250.0
        T2 = 350.0
        iceflag = 3
        with self.assertRaises(ValueError):
            compute_rsat(T, p, T1, T2, iceflag)



class TestMoislif(unittest.TestCase):
    def test_typical_input(self):
        T = 288.15  # K
        qv = 0.014  # kg/kg
        qvv = 0.014  # kg/kg
        qvi = 0.0  # kg/kg
        p0 = 101325  # Pa
        T0 = 288.15  # K
        q0 = 0.01  # kg/kg
        qt = 0.01  # kg/kg
        fracent = 0.001  # m^-1
        prate = 0.0  # m^-1
        T1 = 273.15  # K
        T2 = 233.15  # K

        gamma_m = moislif(T, qv, qvv, qvi, p0, T0, q0, qt, fracent, prate, T1, T2)
        self.assertLess(gamma_m, 0.0)  # lapse rate should be positive

    def test_upperair_input(self):
        T = 260.0  # K
        qv = 0.001  # kg/kg
        qvv = 0.001  # kg/kg
        qvi = 0.001  # kg/kg
        p0 = 70000  # Pa
        T0 = 256.15  # K
        q0 = 0.001  # kg/kg
        qt = 0.01  # kg/kg
        fracent = 0.001  # m^-1
        prate = 0.0  # m^-1
        T1 = 273.15  # K
        T2 = 233.15  # K

        gamma_m = moislif(T, qv, qvv, qvi, p0, T0, q0, qt, fracent, prate, T1, T2)
        self.assertLess(gamma_m, 0.0)  # lapse rate should be positive

    """
    def test_invalid_input(self):
        T = -100.0  # K
        qv = 0.01  # kg/kg
        qvv = 0.01  # kg/kg
        qvi = 0.0  # kg/kg
        p0 = 101325  # Pa
        T0 = 288.15  # K
        q0 = 0.01  # kg/kg
        qt = 0.01  # kg/kg
        fracent = 0.001  # m^-1
        prate = 0.0  # m^-1
        T1 = 273.15  # K
        T2 = 233.15  # K

        with self.assertRaises(ValueError):
            moislif(T, qv, qvv, qvi, p0, T0, q0, qt, fracent, prate, T1, T2)
    """

if __name__ == '__main__':
    unittest.main()
