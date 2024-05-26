#!/usr/bin/python3
import numpy as np
import src.meteolib as meteolib
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# ---------------------------------------------------------------------------------------------------------------------


def create_stuve():
    _, ax = plt.subplots(1, figsize=(15, 11))

    ax.grid(True)

    pmax = 1000
    pmin = 100
    dp = -10
    presvals = np.arange(pmax, pmin+dp, dp)

    # plot moist-adiabats
    for t in np.array([-40, -30, -20, -10, 0, 5, 10, 15, 20, 25, 30, 35, 40]):
        tw = []
        for p in presvals:
            tw.append(meteolib.wetlift(1000., t, p))
        ax.semilogy(tw, presvals, 'b--', lw=0.5, alpha=0.7)

    # plot dry adiabats
    for t in np.array([-80, -60, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230]):
        ax.semilogy(meteolib.thetas(t, presvals), presvals, 'r--', lw=0.5, alpha=0.7)

    # plot lines of mixing ratio
    pmin = 400
    presvals = np.arange(pmax, pmin+dp, dp)

    mrr = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0])
    for i in range(0, np.size(mrr)):
        ax.semilogy(meteolib.temp_at_mixrat(mrr[i], presvals), presvals, 'm--', lw=0.6, alpha=0.8)

    plt.yscale('log')

    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks(np.linspace(100, 1000, 10))
    ax.set_ylim(1050, 400)
    ax.set_xticks(np.arange(-90, 50, 10))
    ax.set_xlim(-90, 50)

    ax.set_title("Stuve Diagramm", fontsize=18)
    plt.xlabel("Temperatur [C]")
    plt.ylabel("Druck [hPa]")
    return ax
