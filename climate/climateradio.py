#!/usr/bin/python3


import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from readhtml import readhtml, col_names

# -----------------------------------------------------------------------------------------------------------------------------
ZEROCNK = 273.15



# -----------------------------------------------------------------------------------------------------------------------------

def plot_tempseries(mean_series, max_series, min_series, std_series, temp, lvl='850', year='2020'):
    _ = plt.figure(figsize=(10, 6))
    ax = plt.subplot(111)
    ax.plot(np.arange(1, 13, 1), mean_series, 'k', lw=2.1, label='Mittel')
    ax.plot(np.arange(1, 13, 1), max_series, 'r--', label='Max')
    ax.plot(np.arange(1, 13, 1), min_series, 'b--', label='Min')
    ax.plot(np.arange(1, 13, 1), temp, 'm', label=year, lw=1.9)

    ax.fill_between(np.arange(1, 13, 1), max_series-0.1, min_series+0.1, where=max_series>= min_series, facecolor='green', alpha=0.4)
    ax.fill_between(np.arange(1, 13, 1), mean_series+std_series, mean_series-std_series, where=mean_series+std_series>= mean_series-std_series, facecolor='green', alpha=0.4)
    if (lvl=='850') or (lvl=='700'):
        ax.axhline(y=273.15, color='k', alpha=0.5)
    elif (lvl=='dt500'):
        ax.axhline(y=25, color='k', alpha=0.5)
    #ax.legend(bbox_to_anchor=(1., 0.75))
    ax.legend(loc=0)
    ax.set_xlim(1, 12)
    ax.set_ylabel('K')
    ax.set_xlabel('Monat')
    if lvl == 'dt500':
        ax.set_title("Lapse Rates 850-500 hPa Monatsmittelwerte von 1961-2014", fontsize=14)
    elif lvl == 'dt300':
        ax.set_title("Lapse Rates 850-300 hPa Monatsmittelwerte von 1961-2014", fontsize=14)
    else:
        ax.set_title(lvl +" hPa Monatsmittelwerte von 1961-2014", fontsize=14)
    plt.savefig(f"hist_{year}_{lvl}.png")
    plt.show()
    plt.close()


# -----------------------------------------------------------------------------------------------------------------------------


def current_year(year=2020, station_id="10393", **kwargs):
    radio_df = pd.DataFrame(columns = col_names)
    today = datetime.date.today()
    month = kwargs.get("month", today.month)
    for month in range(1, month+1):
        for day in range(1, 32):
            for time in [0, 6, 12, 18]:
                datum_liste = (year, month, day, time)
                urlstring = (f"http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&"
                             f"YEAR={year:04d}&MONTH={month:02d}&FROM={day:02d}{time:02d}&TO={day:02d}{time:02d}&STNM={station_id}"
                            )
                radio_df = pd.concat([radio_df, readhtml(urlstring, datum_liste)], ignore_index=True)
    radio_df.to_pickle(f"./climate_{station_id}_{year}.pkl")
    
    temp850 = np.full(12, np.nan)
    temp700 = np.full(12, np.nan)
    temp500 = np.full(12, np.nan)
    temp400 = np.full(12, np.nan)
    temp300 = np.full(12, np.nan)
    temp250 = np.full(12, np.nan)
    temp200 = np.full(12, np.nan)
    temp150 = np.full(12, np.nan)
    temp100 = np.full(12, np.nan)
    dt500 = np.full(12, np.nan)
    dt300 = np.full(12, np.nan)

    for i in range(1, month+1):
        temp850[i-1] = radio_df[(radio_df['pres'] == 850.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        temp700[i-1] = radio_df[(radio_df['pres'] == 700.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        temp500[i-1] = radio_df[(radio_df['pres'] == 500.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        temp400[i-1] = radio_df[(radio_df['pres'] == 400.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        temp300[i-1] = radio_df[(radio_df['pres'] == 300.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        temp250[i-1] = radio_df[(radio_df['pres'] == 250.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        temp200[i-1] = radio_df[(radio_df['pres'] == 200.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        temp150[i-1] = radio_df[(radio_df['pres'] == 150.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        temp100[i-1] = radio_df[(radio_df['pres'] == 100.0) & ((radio_df['month'] == i))]['temp'].mean() + ZEROCNK
        dt500[i-1] = radio_df[(radio_df['pres'] == 850.0) & ((radio_df['month'] == i))]['temp'].mean() - radio_df[(radio_df['pres'] == 500.0) & ((radio_df['month'] == i))]['temp'].mean()
        dt300[i-1] = radio_df[(radio_df['pres'] == 850.0) & ((radio_df['month'] == i))]['temp'].mean() - radio_df[(radio_df['pres'] == 300.0) & ((radio_df['month'] == i))]['temp'].mean()

    datafile = open('10393.csv', 'r').read()
    data = np.array([l.strip() for l in datafile.split('\n')])
    st_id, years, month, hPa850, hPa700, hPa500, hPa400, hPa300, hPa250, hPa200, hPa150, hPa100 = np.loadtxt(data, delimiter=';', unpack=True )
    hPa200[np.where(hPa200==1.0e+20)] = np.nan
    hPa300[np.where(hPa300==1.0e+20)] = np.nan
    hPa500[np.where(hPa500==1.0e+20)] = np.nan
    hPa400[np.where(hPa400==1.0e+20)] = np.nan
    hPa700[np.where(hPa700==1.0e+20)] = np.nan
    hPa850[np.where(hPa850==1.0e+20)] = np.nan
    hPa250[np.where(hPa250==1.0e+20)] = np.nan
    hPa200[np.where(hPa200==1.0e+20)] = np.nan
    hPa150[np.where(hPa150==1.0e+20)] = np.nan
    hPa100[np.where(hPa100==1.0e+20)] = np.nan

    max850, m850, min850, std850 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    max700, m700, min700, std700 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    max500, m500, min500, std500 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    max400, m400, min400, std400 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    max300, m300, min300, std300 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    max250, m250, min250, std250 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    max200, m200, min200, std200 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    max150, m150, min150, std150 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    max100, m100, min100, std100 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    maxdt500, mdt500, mindt500, stddt500 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)
    maxdt300, mdt300, mindt300, stddt300 = np.zeros(12), np.zeros(12), np.zeros(12), np.zeros(12)


    for m in range(1, 13):
        m850[m-1] = np.nanmean(hPa850[np.where(month==m)])
        max850[m-1] = np.nanmax(hPa850[np.where(month==m)])
        min850[m-1] = np.nanmin(hPa850[np.where(month==m)])
        std850[m-1] = np.nanstd(hPa850[np.where(month==m)])
        m250[m-1] = np.nanmean(hPa250[np.where(month==m)])
        max250[m-1] = np.nanmax(hPa250[np.where(month==m)])
        min250[m-1] = np.nanmin(hPa250[np.where(month==m)])
        std250[m-1] = np.nanstd(hPa250[np.where(month==m)])
        m150[m-1] = np.nanmean(hPa150[np.where(month==m)])
        max150[m-1] = np.nanmax(hPa150[np.where(month==m)])
        min150[m-1] = np.nanmin(hPa150[np.where(month==m)])
        std150[m-1] = np.nanstd(hPa150[np.where(month==m)])
        m700[m-1] = np.nanmean(hPa700[np.where(month==m)])
        max700[m-1] = np.nanmax(hPa700[np.where(month==m)])
        min700[m-1] = np.nanmin(hPa700[np.where(month==m)])
        std700[m-1] = np.nanstd(hPa700[np.where(month==m)])
        m500[m-1] = np.nanmean(hPa500[np.where(month==m)])
        max500[m-1] = np.nanmax(hPa500[np.where(month==m)])
        min500[m-1] = np.nanmin(hPa500[np.where(month==m)])
        std500[m-1] = np.nanstd(hPa500[np.where(month==m)])
        m400[m-1] = np.nanmean(hPa400[np.where(month==m)])
        max400[m-1] = np.nanmax(hPa400[np.where(month==m)])
        min400[m-1] = np.nanmin(hPa400[np.where(month==m)])
        std400[m-1] = np.nanstd(hPa400[np.where(month==m)])
        m300[m-1] = np.nanmean(hPa300[np.where(month==m)])
        max300[m-1] = np.nanmax(hPa300[np.where(month==m)])
        min300[m-1] = np.nanmin(hPa300[np.where(month==m)])
        std300[m-1] = np.nanstd(hPa300[np.where(month==m)])
        m200[m-1] = np.nanmean(hPa200[np.where(month==m)])
        max200[m-1] = np.nanmax(hPa200[np.where(month==m)])
        min200[m-1] = np.nanmin(hPa200[np.where(month==m)])
        std200[m-1] = np.nanstd(hPa200[np.where(month==m)])
        m100[m-1] = np.nanmean(hPa100[np.where(month==m)])
        max100[m-1] = np.nanmax(hPa100[np.where(month==m)])
        min100[m-1] = np.nanmin(hPa100[np.where(month==m)])
        std100[m-1] = np.nanstd(hPa100[np.where(month==m)])
        mdt500[m-1] = np.nanmean(hPa850[np.where(month==m)]-hPa500[np.where(month==m)])
        maxdt500[m-1] = np.nanmax(hPa850[np.where(month==m)]-hPa500[np.where(month==m)])
        mindt500[m-1] = np.nanmin(hPa850[np.where(month==m)]-hPa500[np.where(month==m)])
        stddt500[m-1] = np.nanstd(hPa850[np.where(month==m)]-hPa500[np.where(month==m)])
        mdt300[m-1] = np.nanmean(hPa850[np.where(month==m)]-hPa300[np.where(month==m)])
        maxdt300[m-1] = np.nanmax(hPa850[np.where(month==m)]-hPa300[np.where(month==m)])
        mindt300[m-1] = np.nanmin(hPa850[np.where(month==m)]-hPa300[np.where(month==m)])
        stddt300[m-1] = np.nanstd(hPa850[np.where(month==m)]-hPa300[np.where(month==m)])


    plot_tempseries(m850, max850, min850, std850, temp850, lvl='850', year=year)
    plot_tempseries(m700, max700, min700, std700, temp700, lvl='700', year=year)
    plot_tempseries(m500, max500, min500, std500, temp500, lvl='500', year=year)
    plot_tempseries(m400, max400, min400, std400, temp400, lvl='400', year=year)
    plot_tempseries(m250, max250, min250, std250, temp250, lvl='250', year=year)
    plot_tempseries(m200, max200, min200, std200, temp200, lvl='200', year=year)
    plot_tempseries(m150, max150, min150, std150, temp100, lvl='150', year=year)
    plot_tempseries(m100, max100, min100, std100, temp100, lvl='100', year=year)
    plot_tempseries(m300, max300, min300, std300, temp300, lvl='300', year=year)
    plot_tempseries(mdt300, maxdt300, mindt300, stddt300, dt300, lvl='dt300', year=year)
    plot_tempseries(mdt500, maxdt500, mindt500, stddt500, dt500, lvl='dt500', year=year)

def create_boxplot():
    datafile = open('10393.csv', 'r').read()
    data = np.array([l.strip() for l in datafile.split('\n')])
    _, _, _, hPa850, hPa700, hPa500, hPa400, hPa300, hPa250, hPa200, hPa150, hPa100 = np.loadtxt(data, delimiter=';', unpack=True )
 
    # create box plot
    plt.boxplot(np.vstack((hPa850[np.where(hPa850!=1.0e+20)], 
                           hPa700[np.where(hPa700!=1.0e+20)],
                           hPa500[np.where(hPa500!=1.0e+20)],
                           hPa300[np.where(hPa300!=1.0e+20)],)).T, labels=['850', '700', '500', '300'])
    plt.savefig("hist_boxplot.png")
    plt.close()


# -----------------------------------------------------------------------------------------------------------------------------

def main():
    create_boxplot()
    #current_year(year=2024, station_id="10393", month=12)

# -----------------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    main()
