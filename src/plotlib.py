#!/usr/bin/python3

import matplotlib
matplotlib.use('Agg')
# Force matplotlib to not use any Xwindows backend.
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from io import StringIO
from matplotlib.projections import register_projection

import numpy as np


from src.utilitylib import station_number2string
import src.meteolib as meteolib
from src.skewT import SkewXAxes
from src.meteolib import cr
# ---------------------------------------------------------------------------------------------------------------------

def add_adiabatic(ax, pmax=1000, pmin=100, dp=-10):

    presvals = np.arange(pmax, pmin+dp, dp)

    # plot moist-adiabats
    for t in np.array([-40, -30, -20, -10, 0, 5, 10, 15, 20, 25, 30, 35, 40]):
        tw = []
        for p in presvals:
            tw.append(meteolib.wetlift(1000., t + meteolib.cr['ZEROCNK'], p) - meteolib.cr['ZEROCNK'])
        ax.semilogy(tw, presvals, 'b--', lw=0.5, alpha=0.7)

    # plot dry adiabats
    for t in np.array([-80, -60, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230]):
        ax.semilogy(meteolib.thetas(t + meteolib.cr['ZEROCNK'], presvals) - meteolib.cr['ZEROCNK'], presvals, 'r--', lw=0.5, alpha=0.7)

    # plot lines of mixing ratio
    pmin = 400
    presvals = np.arange(pmax, pmin+dp, dp)

    mrr = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0])
    for i in range(0, np.size(mrr)):
        ax.semilogy(meteolib.temp_at_mixrat(mrr[i], presvals), presvals, 'm--', lw=0.6, alpha=0.8)

    return ax

# ---------------------------------------------------------------------------------------------------------------------


def create_stuve():
    fig, ax = plt.subplots(1, figsize=(15, 11))
    plt.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.96)

    ax.grid(True)

    plt.yscale('log')

    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks(np.linspace(100, 1000, 10))
    ax.set_ylim(1050, 400)
    ax.set_xticks(np.arange(-90, 50, 10))
    ax.set_xlim(-90, 50)

    plt.xlabel("Temperatur [C]")
    plt.ylabel("Druck [hPa]")
    return fig, ax

def plot_stuve_cm1(t_env, td_env, pres_env, T_lift, Qv_lift, Qt_lift, T_lif_ecape, Qv_lif_ecape, Qt_lif_ecape):
    _, ax = create_stuve()

    ax = add_adiabatic(ax)

    # plot souding datas
    ax.plot(t_env, pres_env/100, '-b', lw=1.5, label="Temperature")
    ax.plot(td_env, pres_env/100, '-g', lw=1.5, label="Dewpoint")

    # plot lifted parcel
    T_rho = T_lift*(1 + (cr['Rv']/cr['Rd'])*Qv_lift - Qt_lift)
    T_rho = T_rho - cr['ZEROCNK']
    ax.plot(T_rho, pres_env/100, 'r--', linewidth=2, label="CAPE Parcel")

    # plot lifted parcel
    T_rho = T_lif_ecape*(1 + (cr['Rv']/cr['Rd'])*Qv_lif_ecape - Qt_lif_ecape)
    T_rho = T_rho - cr['ZEROCNK']
    ax.plot(T_rho, pres_env/100, 'm--', linewidth=2, label="ECAPE Parcel")

    plt.yscale('log')

    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks(np.linspace(100, 1000, 10))
    ax.set_ylim(1020, 100)
    ax.set_xticks(np.arange(-80, 40, 10))
    ax.set_xlim(-80, 40)

    plt.legend()

    plt.savefig(f"ecape.png")
    plt.close()

def plot_stuve(station_sounding_obj, number_str):
    _, ax = create_stuve()

    ax = add_adiabatic(ax)

    # plot souding datas
    t_env = station_sounding_obj.get_sounding_temp()
    tvir_env = station_sounding_obj.get_sounding_tv()
    pres_env = station_sounding_obj.get_sounding_pres()
    ax.plot(station_sounding_obj.get_sounding_wetbulb(), pres_env, '-c', lw=1.1)
    ax.plot(tvir_env, pres_env, '-b', lw=1.1)
    ax.plot(t_env, pres_env, '-b', lw=1.5)
    ax.plot(station_sounding_obj.get_sounding_dew(), pres_env, '-g', lw=1.5)

    t_parcel = station_sounding_obj.SB_parcel - meteolib.cr['ZEROCNK']

    ax.plot(t_parcel, pres_env, '-r', lw=1.8, label="SB CAPE Parcel")

    if station_sounding_obj.SB_CAPE > 10:
        ax.fill_betweenx(pres_env, tvir_env, t_parcel, where=t_parcel >= tvir_env, facecolor='red', alpha=0.4)

    t_parcel = station_sounding_obj.ECAPE_parcel - meteolib.cr['ZEROCNK']

    ax.plot(t_parcel, pres_env, 'm--', linewidth=2, label="ECAPE Parcel")
    if station_sounding_obj.ECAPE > 10:
        ax.fill_betweenx(pres_env, tvir_env, t_parcel, where=t_parcel >= tvir_env, facecolor='red', alpha=0.2)
    
    plt.yscale('log')

    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks(np.linspace(100, 1000, 10))
    ax.set_ylim(1050, 100)
    ax.set_xticks(np.arange(-90, 50, 10))
    ax.set_xlim(-90, 50)

    #ax.axhline(station_sounding_obj.efibot)
    #ax.axhline(station_sounding_obj.efitop)

    # wind bars
    anzahl_bar = 45
    anzahl = (np.where(station_sounding_obj.pres_env <= 100))[0][0]
    iter_bar = int(anzahl/anzahl_bar)
    x_value = np.zeros(np.size(station_sounding_obj.pres_env[:anzahl:iter_bar])) + 45
    ax.barbs(x_value, station_sounding_obj.pres_env[:anzahl:iter_bar],
             station_sounding_obj.u_env[:anzahl:iter_bar], station_sounding_obj.v_env[:anzahl:iter_bar],
             station_sounding_obj.wind_speed[:anzahl:iter_bar],
             length=8, pivot='middle', barb_increments=dict(half=2.5, full=5, flag=25))

    ax.set_title(f"{station_number2string(number_str)}", fontsize=18)
    plt.legend()
    plt.savefig(f"{station_number2string(number_str)}_stuve.png")
    plt.close()

# ---------------------------------------------------------------------------------------------------------------------
# Create a SkewT Plots


def plot_skewT(station_sounding_obj, number_str):

    register_projection(SkewXAxes)

    # fig = plt.figure(figsize=(6.5875, 6.2125))
    fig = plt.figure(figsize=(9, 6.2125))
    ax = fig.add_subplot(111, projection='skewx')
    plt.subplots_adjust(left=0.11, right=0.7, bottom=0.07, top=0.91)
    # Grid
    ax.grid(True)

    # add lines
    #ax = add_adiabatic(ax)

    # lines
    pmax = 1000
    pmin = 100
    dp = -10
    presvals = np.arange(pmax, pmin+dp, dp)

    # plot moist-adiabats
    for t in np.array([-40, -30, -20, -10, 0, 5, 10, 15, 20, 25, 30, 35, 40]):
        tw = []
        for p in presvals:
            tw.append(meteolib.wetlift(1000., t + meteolib.cr['ZEROCNK'], p) - meteolib.cr['ZEROCNK'])
        ax.semilogy(tw, presvals, 'b--', lw=0.7, alpha=0.7)

    # plot dry adiabats
    for t in np.array([-80, -60, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 70, 90, 110, 130, 150, 170, 190, 210, 230]):
        ax.semilogy(meteolib.thetas(t + meteolib.cr['ZEROCNK'], presvals) - meteolib.cr['ZEROCNK'], presvals, 'r--', lw=0.7, alpha=0.7)

    # lines of mixing ratio
    pmin = 700
    presvals = np.arange(pmax, pmin+dp, dp)

    mrr = np.append(station_sounding_obj.mr_env[0]*1000, np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]))
    for i in range(0, np.size(mrr)):
        ax.semilogy(meteolib.temp_at_mixrat(mrr[i], presvals), presvals, 'm--', lw=0.6, alpha=0.8)

    # plot tile/location
    plt.annotate(f"{station_number2string(number_str)}", xy=(0.25, 1.04), xycoords='axes fraction',fontsize=15)
    plt.xlabel("Temperature [C]")
    plt.ylabel("Pressure [hPa]")

    ax.axvline(0, color='b', linestyle='--' , lw=0.9)
    ax.axvline(-20, color='b', linestyle='--', lw=0.9)

    t_env = station_sounding_obj.get_sounding_temp()
    tvir_env = station_sounding_obj.get_sounding_tv()
    pres_env = station_sounding_obj.get_sounding_pres()
    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    ax.semilogy(station_sounding_obj.get_sounding_wetbulb(), pres_env, 'c-', lw=0.8)  # Plot the wetbulb profile
    ax.semilogy(tvir_env,pres_env,'b-'  , lw=1.1, alpha=0.9)  # Plot the virtuell temperature
    ax.semilogy(t_env, pres_env, 'b-' , lw=1.8 )  # Plot the temperature profile
    ax.semilogy(station_sounding_obj.get_sounding_dew(), pres_env, 'g-', lw=1.8)  # plot the dewpoint profile

    t_parcel = station_sounding_obj.SB_parcel - meteolib.cr['ZEROCNK']
    ax.semilogy(t_parcel, pres_env, '-r', lw=1.4)
    ax.fill_betweenx(pres_env, tvir_env, t_parcel,  where=t_parcel > tvir_env, facecolor='red', interpolate=True, alpha=0.4)

    t_parcel = station_sounding_obj.ECAPE_parcel - meteolib.cr['ZEROCNK']

    ax.semilogy(t_parcel, pres_env, 'm--', linewidth=2, label="ECAPE Parcel")
    ax.fill_betweenx(pres_env, tvir_env, t_parcel, where=t_parcel >= tvir_env, facecolor='red', alpha=0.2)

    # Disables the log-formatting that comes with semilogy
    ax.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_yticks(np.linspace(100, 1000, 10))
    ax.set_ylim(1050, 100)

    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    #ax.set_xticks(np.arange(-30, 50, 10))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[60] = ' '
    labels[59] = ' '
    labels[:53] = ' '
    ax.set_xticklabels(labels)
    ax.set_xlim(-35, 40)

    # Draw the hodograph on the Skew-T.
    ax2 = plt.axes([.725, .40, .25, .50])

    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    clim=20
    if clim > 20:
        for i in range(10, clim+20, 10):
            # Draw the range rings around the hodograph.
            circle = plt.Circle((0, 0), i, color='k', alpha=.3, fill=False)
            ax2.add_artist(circle)
    else:
        for i in range(5, 40, 5):
            # Draw the range rings around the hodograph.
            circle = plt.Circle((0, 0), i, color='k', alpha=.3, fill=False)
            ax2.add_artist(circle)
    height = station_sounding_obj.get_sounding_height()
    i = 0
    while height[i] < 10000 :
            i += 1
    points = np.array([station_sounding_obj.u_env[:i+1], station_sounding_obj.v_env[:i+1]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # Create a continuous norm to map from data points to colors
    # cmap = ListedColormap(['r', 'g', 'b'])
    norm = plt.Normalize(0, 10000)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red", "purple", "navy", "forestgreen", "darkgreen"])
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    # Set the values used for colormapping

    lc.set_array(height)
    lc.set_linewidth(2)
    line = ax2.add_collection(lc)
    cbar = fig.colorbar(line, ax=ax2, orientation='horizontal', ticks=np.arange(0, 10000, 1000))
    clevs = np.arange(0,10,1)
    cbar.ax.set_xticklabels(labels=clevs[::1])
    mw6 = station_sounding_obj.get_mw_0_6km()

    rstu, rstv, lstu, lstv = station_sounding_obj.get_sorm_motion()
    ax2.plot(rstu, rstv, 'ro', alpha=0.9, zorder=2)  # Plot Bunker's Storm motion right mover as a red dot
    ax2.plot(lstu, lstv, 'bo', alpha=0.9, zorder=2)  # Plot Bunker's Storm motion left mover as a blue dot
    ax2.plot(mw6[0], mw6[1], 'o', color = 'grey', alpha=0.7, zorder=2)  # 0-6km Mean weight wind

    ax2.axhline(y=0, color='k', zorder=1, lw=1.1)
    ax2.axvline(x=0, color='k', zorder=1, lw=1.1)

    ax2.set_xlim(-clim, clim)
    ax2.set_ylim(-clim, clim)
    if clim == 40:
        ax2.text(-5.2, -39.7, "40 m/s")
        ax2.text(-5.2, -29.7, "30")
        ax2.text(-5.2, -19.7, "20")
        #ax2.text(-4.82, -9.7, "10 m/s")
    elif clim == 20:
        #ax2.text(-4.43, -29.7, "30 m/s")
        ax2.text(-3.6, -19.7, "20 m/s")
        ax2.text(-3.6, -9.7, "10 m/s")
    elif clim == 50:
        ax2.text(-6.2, -49.7, "50 m/s")
        ax2.text(-6.2, -39.7, "40")
        ax2.text(-6.2, -29.7, "30")
    else:
        ax2.text(-4.48, -29.7, "30 m/s")
        ax2.text(-4.48, -19.7, "20 m/s")
        ax2.text(-4.48,-9.7,"10 m/s")
    fig.text(0.850, 0.39, "km")
    x_value = 0.725
    mw_dir, mw_speed = meteolib.uv2spddir(mw6[0], mw6[1])
    fig.text(x_value, 0.33, r"Mean Wind (0-6km) : %3.0f$^{\circ}$/%.0f m/s" % (mw_dir, mw_speed))
    mw8 = station_sounding_obj.get_mw_0_8km()
    mw_dir, mw_speed = meteolib.uv2spddir(mw8[0], mw8[1])
    fig.text(x_value, 0.30, r"Mean Wind (0-8km) : %3.0f$^{\circ}$/%.0f m/s" % (mw_dir, mw_speed))

    #fig.text(x_value, 0.30, "ML CAPE : %.1f J/kg" % (2843.2))
    fig.text(x_value, 0.27, "SB CAPE : %.1f J/kg" % (station_sounding_obj.get_sb_cape()))
    fig.text(x_value, 0.24, "ECAPE   : %.1f J/kg" % (station_sounding_obj.ECAPE))

    #fig.text(x_value, 0.24, "BRN ML : %.1f" % (brn))
    wmaxshr = station_sounding_obj.get_wmaxshear()
    if (wmaxshr >= 1000.0):
        fig.text(x_value,0.21, r"WMAXSHEAR ML : %.1f $m^2/s^2$" % (wmaxshr), color='r', weight='bold')
    elif (wmaxshr < 1000.0) and (wmaxshr >= 400.0):
        fig.text(x_value,0.21, r"WMAXSHEAR ML : %.1f $m^2/s^2$" % (wmaxshr), color='orange', weight='bold')
    else:
        fig.text(x_value,0.21, r"WMAXSHEAR ML : %.1f $m^2/s^2$" % (wmaxshr))
    """
    if pw >= 40.0:
        fig.text(x_value, 0.15,r"PW : %.1f $kg/m^2$" % (pw), color='r', weight='bold')
    elif (pw < 40.0) and (pw > 30.0):
        fig.text(x_value, 0.15,r"PW : %.1f $kg/m^2$" % (pw), color='orange', weight='bold')
    else:
        fig.text(x_value, 0.15,r"PW : %.1f $kg/m^2$" % (pw))
    """
    #fig.text(x_value, 0.18, "Specific humidity (100mb) : %.1f g/kg" % (q_mean))
    #fig.text(x_value, 0.15, "Lapse Rate (850-500hPa) : %.1f K/km" %(lapse))
    try:
        plt.savefig(f"{station_number2string(number_str)}_skewT.png")
    except ValueError:
        pass
    plt.close()
