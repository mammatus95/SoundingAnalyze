#!/usr/bin/python3
import sys

# project moduls
from src.sounding import sounding
from src.plotlib import plot_stuve
# ---------------------------------------------------------------------------------------------------------------------

def main():
    urlstring = (f"http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime={sys.argv[1]}%20{sys.argv[2]}:00:00"
                 f"&id={sys.argv[3]}&type=TEXT:LIST")

    station_sounding_obj = sounding(urlstring)
    print(station_sounding_obj.get_meanwind())
    print("Surface CAPE", station_sounding_obj.get_sb_cape())
    plot_stuve(station_sounding_obj, sys.argv[3])

if __name__ == "__main__":
    main()
