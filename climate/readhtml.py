#!/usr/bin/python3

import time
import urllib.request
import pandas as pd


# -----------------------------------------------------------------------------------------------------------------------------


col_names = ['id', 'pres', 'year', 'month', 'day', 'time', 'temp', 'height', 'dew', 'mix', 'wdir', 'wspd', 'pot']

def KT2MS(v):
    return v*0.514444

# -----------------------------------------------------------------------------------------------------------------------------
# read wyoming data


def data2df(data, station_id, datum_liste):
    """
    -----------------------------------------------------------------------------
    PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
    hPa     m      C      C      %    g/kg    deg   knot     K      K      K
    -----------------------------------------------------------------------------
    """
    row = [[int(station_id), float(data[0]), datum_liste[0], datum_liste[1], datum_liste[2], datum_liste[3], 
            float(data[2]), float(data[1]), float(data[3]), float(data[5]), int(data[6]), 
            KT2MS(float(data[7])), float(data[8])]]
    return row


def readhtml(urlstring, datum_liste, station_id=10393, mode=1):

    radio_df = pd.DataFrame(columns = col_names)
    pres_list = ['1000.0', '850.0', '700.0', '600.0', '500.0', '400.0', '300.0', '250.0', '200.0', '150.0', '100.0']

    try:
        fp = urllib.request.urlopen(urlstring)
    except urllib.error.HTTPError as e:
        if (e.code == 503):
            print(f"Error-Code: {e.code} {e.read()}")
            time.sleep(10)
            fp = urllib.request.urlopen(urlstring)
        else:
            print(str(e.code), " ", e.read())
            print(urlstring)
            return radio_df
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()

    seq = mystr.split("</PRE>")
    if mode == 1:
        mystr = seq[0]
        liste = seq[0].splitlines(True)
        for i in range(10, len(liste)-2):
            data=[]
            if len(liste[i].split()) > 10:
                data = liste[i].split()
                if data[0] in pres_list:
                    df = pd.DataFrame(data2df(data, station_id, datum_liste), 
                                      columns = col_names)
                    radio_df = pd.concat([radio_df, df], ignore_index=True)
    else:
        liste = seq[2].splitlines(True)
        for i in range(6, len(liste)-1):
            if len(liste[i].split()) > 6:
                print(liste[i])
    return radio_df

# -----------------------------------------------------------------------------------------------------------------------------


def test_read_wyoming(year=1992, month=4, day=20, time=12, station_id="10393"):

    radio_df = pd.DataFrame(columns = col_names)
    datum_liste = (year, month, day, time)
    urlstring = "http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%%3ALIST&YEAR=%.4d&MONTH=%.2d&FROM=%.2d%.2d&TO=%.2d%.2d&STNM=%s" % (year, month, day, time, day, time, station_id)
    radio_df = radio_df.append(readhtml(urlstring, datum_liste), ignore_index=True)
    print(radio_df)
    # check instance of pandas dataframe
    assert isinstance(radio_df, pd.DataFrame) 
    assert radio_df.shape[1] == len(col_names)


# -----------------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    test_read_wyoming()
