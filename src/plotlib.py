#!/usr/bin/python3

import matplotlib
matplotlib.use('Agg')
# Force matplotlib to not use any Xwindows backend.
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from io import StringIO

import numpy as np
import numpy.ma as ma
import sys
import time
from copy import deepcopy

from src.utilitylib import station_number2strig
import src.meteolib as meteolib
from src.skewT import *


def plot_stuve(station_sounding_obj, number_str):
    b, ax = plt.subplots(1, figsize=(15, 11))
    plt.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.96)
    ax.grid(True)

    pmax = 1000
    pmin = 100
    dp = -10
    presvals = np.arange(pmax, pmin+dp, dp)

    # plot the moist-adiabats
    for t in np.array([-40, -30, -20, -10, 0, 5, 10, 15, 20, 25, 30, 35, 40]):
        tw = []
        for p in presvals:
            tw.append(meteolib.wetlift(1000., t + meteolib.cr['ZEROCNK'], p))
        ax.semilogy(tw, presvals, 'b--', lw=0.7, alpha=0.7)

    # plot the dry adiabats
    for t in np.array([-80, -60, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230]):
        ax.semilogy(meteolib.thetas(t+meteolib.cr['ZEROCNK'], presvals)-meteolib.cr['ZEROCNK'], presvals, 'r--', lw=0.7, alpha=0.7)

    # linie gleichensaettiungsmischungsverhaeltnis
    pmin = 700
    presvals = np.arange(pmax, pmin+dp, dp)

    mrr = np.append(station_sounding_obj.mr_env[0]*1000, np.array([0.1, 0.5, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0]))
    for i in range(0, np.size(mrr)):
        ax.semilogy(meteolib.temp_at_mixrat(mrr[i], presvals), presvals, 'm--', lw=0.6, alpha=0.8)

    # plot souding datas
    print(station_sounding_obj.wetbulb_env - meteolib.cr['ZEROCNK'])
    ax.plot(station_sounding_obj.wetbulb_env - meteolib.cr['ZEROCNK'], station_sounding_obj.pres_env, '-c', lw=1.1)
    ax.plot(station_sounding_obj.tvir_env - meteolib.cr['ZEROCNK'], station_sounding_obj.pres_env, '-b', lw=1.1)
    ax.plot(station_sounding_obj.t_env - meteolib.cr['ZEROCNK'], station_sounding_obj.pres_env, '-b', station_sounding_obj.dew_env - meteolib.cr['ZEROCNK'], station_sounding_obj.pres_env, '-g', lw=1.5)
    #ax.plot(station_sounding_obj.mu_parcel, station_sounding_obj.press_parcel, '-r', lw=1.8)
    #ax.plot(mu2_parcel, press_parcel2, '-m', lw=1.9)
    """
    if float(sys.argv[3]) > 10.0:
        mu = np.zeros(station_sounding_obj.pres_env.size, dtype=np.float64)
        for i in range(int(sys.argv[2]), station_sounding_obj.pres_env.size):
            #if pres_env[i] >= lclp:
            #    mu[i] = thetas(theta(pres_env[int(sys.argv[2])], t_temp[int(sys.argv[2])], p2=1000.), pres_env[i])
            #else:
            mu[i] = meteolib.wetlift(station_sounding_obj.lclp, station_sounding_obj.lclt, station_sounding_obj.pres_env[i])

        ax.fill_betweenx(station_sounding_obj.pres_env, station_sounding_obj.t_temp, mu, where=mu >= station_sounding_obj.t_temp, facecolor='red', alpha=0.4)
    """

    plt.yscale('log')

    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks(np.linspace(100, 1000, 10))
    ax.set_ylim(1050, 100)
    ax.set_xticks(np.arange(-90, 50, 10))
    ax.set_xlim(-90, 50)

    #plt.annotate(string, xy=(float(sys.argv[5])+20, float(sys.argv[4])), xytext=(float(sys.argv[5]) - 85, float(sys.argv[4])+10), arrowprops=dict(arrowstyle="-", facecolor='black'))
    #ax.axhline(float(sys.argv[10])) # efi bot
    #ax.axhline(float(sys.argv[11])) # efi top

    """
    # wind bars
    anzahl_bar = 45
    anzahl = (np.where(station_sounding_obj.pres_env==100))[0][0]
    iter_bar = int(anzahl/anzahl_bar)
    x = np.zeros(np.size(station_sounding_obj.pres_env[:anzahl:iter_bar])) + 45
    ax.barbs(x, station_sounding_obj.pres_env[:anzahl:iter_bar], 
             station_sounding_obj.u_env[:anzahl:iter_bar], station_sounding_obj.v_env[:anzahl:iter_bar], 
             station_sounding_obj.wind_speed, length=8, pivot='middle', barb_increments=dict(half=2.5, full=5, flag=25))
    """
    #plt.text(-90, 99, wenn, fontsize=11)
    #ax.set_title(title, fontsize=18)
    plt.xlabel("Temperatur [C]")
    plt.ylabel("Druck [hPa]")
    plt.savefig(f"{station_number2strig(number_str)}_thermo.png")
    plt.close()

