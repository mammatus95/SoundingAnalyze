#!/usr/bin/python3
import unittest
import numpy as np
# project modul
from src.meteolib import cr, temp_at_mixrat, thalvl, theta, thetas, uvwind, uv2spddir

# ----------------------------------------------------------------------------------------------------------------------------


class TestTempAtMixrat(unittest.TestCase):
    def test_valid_inputs(self):
        w = np.array([1, 2, 3])  # mixing ratio (g/kg)
        p = np.array([1000, 900, 800])  # pressure (hPa)
        expected_result = np.array([-17.0913, -9.89, -6.22])  # expected temperature (C)
        result = temp_at_mixrat(w, p)
        np.testing.assert_allclose(result, expected_result, atol=0.01)

    def test_edge_cases(self):
        w = np.array([30, 1e-6])  # edge cases for mixing ratio
        p = np.array([1000, 1000])  # constant pressure
        result = temp_at_mixrat(w, p)
        #self.assertTrue(np.isnan(result[0]))  # expect NaN for zero mixing ratio
        #self.assertGreater(result[1], -100)  # expect finite temperature for small mixing ratio


class TestThalvlFunction(unittest.TestCase):

    def test_valid_input(self):
        theta = np.array([300, 310, 320])
        t = np.array([280, 290, 300])
        result = thalvl(theta, t)
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, theta.shape)

    def test_edge_cases(self):
        theta = np.array([300, 300, 10])
        t = np.array([300, 10, 300])
        result = thalvl(theta, t)
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, theta.shape)


class TestThetaFunction(unittest.TestCase):

    def test_valid_inputs_default_p0(self):
        pres = np.array([1000, 900, 800])
        temp = np.array([300, 310, 320])
        expected_result = np.array([300, 319.47381, 341.066097])
        result = theta(pres, temp)
        np.testing.assert_allclose(result, expected_result, atol=0.01)

    def test_multiple_input_values(self):
        pres = np.array([1000, 900, 800, 700, 600])
        temp = np.array([300, 310, 320, 330, 340])
        expected_result = np.array([300, 319.47, 341.07, 365.40, 393.43])
        result = theta(pres, temp)
        np.testing.assert_allclose(result, expected_result, atol=0.01)


class TestThetasFunction(unittest.TestCase):

    def test_default_p0(self):
        theta = 300
        presvals = 1000
        expected_result = theta
        self.assertAlmostEqual(thetas(theta, presvals), expected_result)

    def test_custom_p0(self):
        theta = 300
        presvals = 1000
        p0 = 500
        expected_result = theta * np.power(presvals/p0, cr['ROCP'])
        self.assertAlmostEqual(thetas(theta, presvals, p0), expected_result)


class TestUVWind(unittest.TestCase):
    def setUp(self):

        # absolute allowed difference esimate and expected value
        self.delta = 0.01 # m/s

    def test_zero_wind_speed(self):
        wind_speed = 0
        wind_dir = 0
        expected_u = 0
        expected_v = 0
        u, v = uvwind(wind_dir, wind_speed)
        self.assertAlmostEqual(u, expected_u, delta=self.delta)
        self.assertAlmostEqual(v, expected_v, delta=self.delta)


    def test_negative_wind_speed(self):
        wind_speed = -10
        wind_dir = 90
        expected_u = -10
        expected_v = 0
        u, v = uvwind(wind_dir, wind_speed)
        self.assertAlmostEqual(u, expected_u, delta=self.delta)
        self.assertAlmostEqual(v, expected_v, delta=self.delta)


    def test_east_wind(self):
        wind_speed = 10
        wind_dir = 90
        expected_u = -10
        expected_v = 0
        u, v = uvwind(wind_dir, wind_speed)
        self.assertAlmostEqual(u, expected_u, delta=self.delta)
        self.assertAlmostEqual(v, expected_v, delta=self.delta)


    def test_west_wind(self):
        wind_speed = 10
        wind_dir = 270
        expected_u = 10
        expected_v = 0
        u, v = uvwind(wind_dir, wind_speed)
        self.assertAlmostEqual(u, expected_u, delta=self.delta)
        self.assertAlmostEqual(v, expected_v, delta=self.delta)


    def test_south_wind(self):
        wind_speed = 10
        wind_dir = 180
        expected_u = 0
        expected_v = 10
        u, v = uvwind(wind_dir, wind_speed)
        self.assertAlmostEqual(u, expected_u, delta=self.delta)
        self.assertAlmostEqual(v, expected_v, delta=self.delta)


    def test_north_wind(self):
        wind_speed = 10
        wind_dir = 0
        expected_u = 0
        expected_v = -10
        u, v = uvwind(wind_dir, wind_speed)
        self.assertAlmostEqual(u, expected_u, delta=self.delta)
        self.assertAlmostEqual(v, expected_v, delta=self.delta)

        wind_dir = 360
        u, v = uvwind(wind_dir, wind_speed)
        self.assertAlmostEqual(u, expected_u, delta=self.delta)
        self.assertAlmostEqual(v, expected_v, delta=self.delta)


class TestUV2SPDDIR(unittest.TestCase):

    def setUp(self):

        # absolute allowed difference esimate and expected value
        self.delta = 0.1 # degree

    def test_scalar_inputs(self):
        u = 1.0
        v = 2.0
        direction, speed = uv2spddir(u, v)
        self.assertIsInstance(direction, float)
        self.assertIsInstance(speed, float)

    def test_numpy_array_inputs(self):
        u = np.array([1.0, 2.0, 3.0])
        v = np.array([4.0, 5.0, 6.0])
        direction, speed = uv2spddir(u, v)
        self.assertIsInstance(direction, np.ndarray)
        self.assertIsInstance(speed, np.ndarray)

    def test_zero_inputs(self):
        u = 0.0
        v = 0.0
        direction, speed = uv2spddir(u, v)
        self.assertTrue(np.isnan(direction))
        self.assertEqual(speed, 0.0)

    def test_east_wind(self):
        u = -10.0
        v = 0.0
        expected_wdir = 90.0
        excepted_speed = 10.0
        direction, speed = uv2spddir(u, v)
        self.assertAlmostEqual(direction, expected_wdir, delta=self.delta)
        self.assertEqual(speed, excepted_speed)

    def test_west_wind(self):
        u = 10.0
        v = 0.0
        expected_wdir = 270.0
        excepted_speed = 10.0
        direction, speed = uv2spddir(u, v)
        self.assertAlmostEqual(direction, expected_wdir, delta=self.delta)
        self.assertEqual(speed, excepted_speed)
    
    def test_south_wind(self):
        u = 0.0
        v = 10.0
        expected_wdir = 180.0
        excepted_speed = 10.0
        direction, speed = uv2spddir(u, v)
        self.assertAlmostEqual(direction, expected_wdir, delta=self.delta)
        self.assertEqual(speed, excepted_speed)

    def test_north_wind(self):
        u = 0.0
        v = -10.0
        expected_wdir = 0.0
        excepted_speed = 10.0
        direction, speed = uv2spddir(u, v)
        self.assertAlmostEqual(direction, expected_wdir, delta=self.delta)
        self.assertEqual(speed, excepted_speed)

    def test_southwest_wind(self):
        u = 3.0
        v = 4.0
        expected_wdir = 216.87
        excepted_speed = 5.0
        direction, speed = uv2spddir(u, v)
        self.assertAlmostEqual(direction, expected_wdir, delta=self.delta)
        self.assertEqual(speed, excepted_speed)

    def test_northwest_wind(self):
        u = 5.0
        v = -3.0
        expected_wdir = 300.96
        direction, _ = uv2spddir(u, v)
        self.assertAlmostEqual(direction, expected_wdir, delta=self.delta)

# ----------------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    unittest.main()
