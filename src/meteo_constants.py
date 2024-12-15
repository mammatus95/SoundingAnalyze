#!/usr/bin/python3
import math
import numpy as np

# -----------------------------------------------------------------------------------------------------------------------------
# meteorology constants

# -----------------------------------------------------------------------------------------------------------------------------
# basic constants

basic_const = {
    'ZEROCNK': 273.15,                       # K        Zero Celsius in Kelvins
    'p0': 100000.0,                          # Pa       reference pressure
    'T0': 273.15,                            # K        reference temperature
    'e0': 611,                               # Pa       reference water vapor (p0, T0)
    'Pas2hPah': 36,                          # Pa/s to hPa/h
    'kn2ms': 0.514,                          # 1 knot = 0.514 m/s
    'ms2kn': 1.94384,                        # 1 m/s  = 1.94384
    'SeaMile2m': 1852.,                      # Sea Miles to m
    'StatuteMile2m': 1609.,                  # Statue Mile to m
    'Re': 6371008.767,                       # m        earth radius
    # 'omega': (2*math.pi)/(23*3600+56*60+4.1),   # rad/s    angular velocity of the earth
    'omega': 7.29211505392569e-05,           # rad/s    angular velocity of the earth
}

# -----------------------------------------------------------------------------------------------------------------------------
# math constants

math_const = {
    'Pi': math.pi,                           # Pi
    'e': 2.71828,                            # Euler Zahl
}

# -----------------------------------------------------------------------------------------------------------------------------
# usage in Peters .et al 2022, 2023

peters_const = {
    'G': 9.81,                               # m/(s*s)  Gravity
    'Rd': 287.04,                            # J/kg/K   dry air gas constant
    'Rv': 461.5,                             # J/kg/K   water vapor gas constant
    'cpd': 1005.7,                           # J/kg/K   isobaric dry air specific heat capacit
    'cpv': 1870.0,                           # J/kg/K   water vapor specific heat capacity (isobaric),
    'cpl': 4190.0,                           # J/kg/K   liquid water specific heat capacity
    'cpi': 2106.0,                           # J/kg/K   ice specific heat capacity
    'Lref': 2.5*1e6,                         # J/kg     latent heat of condensation at reference conditions
    'xlv': 2501000,                          # J/kg     reference latent heat of vaporization at the triple point temperature
    'xls': 2834000,                          # J/kg     reference latent heat of sublimation at the triple point temperature
    'eref': 611.25,                          # Pa       partial pressure of water vapor a triple point temp
    'ttrip': basic_const['T0'],              # K        triple point temp
    'ksq': 0.18,                             # von Karman constant
    'Pr': 1/3,                               # Prandtl number
    'Lm': 250,                               # horizontal mixing Length
    'c1': 0.0498646455,
    'c2': 2.4082965,
    'c3': 7.07475,
    'c4': 38.9114,
    'c5': 0.0915,
    'c6': 1.2035
}

# -----------------------------------------------------------------------------------------------------------------------------
# usage in Bolton 1980

bolton_const = {
    'G': 9.80665,                            # m/(s*s)  Gravity
    'Rd': 287.04,                            # J/kg/K   dry air gas constant
    'Rv': 461.5,                             # J/kg/K   water vapor gas constant
    'cpd': 1005.7,                           # J/kg/K   isobaric dry air specific heat capacit
    'cpv': 1870.0,                           # J/kg/K   water vapor specific heat capacity (isobaric),
    'cpl': 4190.0,                           # J/kg/K   liquid water specific heat capacity
    'cpi': 2106.0,                           # J/kg/K   ice specific heat capacity
    'Lref': 2.5*1e6,                         # J/kg     latent heat of condensation at reference conditions
    'xlv': 2501000,                          # J/kg     reference latent heat of vaporization at the triple point temperature
    'xls': 2834000,                          # J/kg     reference latent heat of sublimation at the triple point temperature
    'eref': 611.25,                          # Pa       partial pressure of water vapor a triple point temp
    'ttrip': basic_const['T0'],              # K        triple point temp
    'ksq': 0.18,                             # von Karman constant
    'Pr': 1/3,                               # Prandtl number
    'Lm': 250,                               # horizontal mixing Length
    'c1': 0.0498646455,
    'c2': 2.4082965,
    'c3': 7.07475,
    'c4': 38.9114,
    'c5': 0.0915,
    'c6': 1.2035
}

# -----------------------------------------------------------------------------------------------------------------------------
# usage in sharppy/nsharp

sharp_const = {
    'G': 9.80665,                            # m/(s*s)  Gravity
    'Rd': 287.04,                            # J/kg/K   dry air gas constant
    'Rv': 461.5,                             # J/kg/K   water vapor gas constant
    'cpd': 1005.7,                           # J/kg/K   isobaric dry air specific heat capacit
    'cpv': 1870.0,                           # J/kg/K   water vapor specific heat capacity (isobaric),
    'cpl': 4190.0,                           # J/kg/K   liquid water specific heat capacity
    'cpi': 2106.0,                           # J/kg/K   ice specific heat capacity
    'epsilon': 0.62197,                      # dimless  molecular weight ratio (Rd/Rv)
    'kappa': 0.28571426,                     # dimless  R/Cp
    'Lref': 2.5*1e6,                         # J/kg     latent heat of condensation at reference conditions
    'xlv': 2501000,                          # J/kg     reference latent heat of vaporization at the triple point temperature
    'xls': 2834000,                          # J/kg     reference latent heat of sublimation at the triple point temperature
    'eref': 611.25,                          # Pa       partial pressure of water vapor a triple point temp
    'ttrip': basic_const['T0'],              # K        triple point temp
    'ksq': 0.18,                             # von Karman constant
    'Pr': 1/3,                               # Prandtl number
    'Lm': 250,                               # horizontal mixing Length
    'c1': 0.0498646455,
    'c2': 2.4082965,
    'c3': 7.07475,
    'c4': 38.9114,
    'c5': 0.0915,
    'c6': 1.2035
}
# -----------------------------------------------------------------------------------------------------------------------------
# put them together

constants = {**basic_const, **math_const, **peters_const, 
             'epsilon': peters_const['Rd']/peters_const['Rv'],
             'kappa': peters_const['Rd']/peters_const['cpd']
             }

constants_bolton = {**basic_const, **math_const, **bolton_const, 
                    'epsilon': bolton_const['Rd']/bolton_const['Rv'],
                    'kappa': bolton_const['Rd']/bolton_const['cpd']
                    }

constants_sharp = {**basic_const, **math_const, **sharp_const}