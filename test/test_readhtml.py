#!/usr/bin/python3
import unittest
import numpy as np
# project modul
from src.readhtml_bufr import readhtml2numpy


class Test_READHTML(unittest.TestCase):
    def setUp(self):
        self.urlstring = "http://weather.uwyo.edu/cgi-bin/bufrraob.py?src=bufr&datetime=2018-12-07%2000:00:00&id=10393&type=TEXT:LIST"
        self.expectedarraysize = 5198

    def test_readhtml2numpy(self):
        self.pres, self.height, self.temp, self.dewpoint, self.mixing_ratio, self.winddir, self.wind_speed = readhtml2numpy(self.urlstring)
        self.assertIsInstance(self.pres, np.ndarray)
        self.assertIsInstance(self.height, np.ndarray)
        self.assertIsInstance(self.temp, np.ndarray)
        self.assertIsInstance(self.dewpoint, np.ndarray)
        self.assertIsInstance(self.mixing_ratio, np.ndarray)
        self.assertIsInstance(self.winddir, np.ndarray)
        self.assertIsInstance(self.wind_speed, np.ndarray)

    def test_pressure(self):
        self.pres, _, _, _, _, _, _ = readhtml2numpy(self.urlstring)
        self.assertGreaterEqual(self.pres.size, self.expectedarraysize, "Pressure Array Size are to small")
        self.assertGreaterEqual(self.pres.max(), 1000, "Surface Pressure is to small")


if __name__ == '__main__':
    unittest.main()
