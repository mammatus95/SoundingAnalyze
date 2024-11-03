#!/usr/bin/python3
import sys

# project moduls
from src.sounding import sounding
from src.plotlib import plot_stuve, plot_skewT
# ---------------------------------------------------------------------------------------------------------------------

def main():
    name = "sounding_base"
    # name = "sounding_noshear"
    sound_filename = f'src/example/{name}.txt'


    station_sounding_obj = sounding(sound_filename)
    print(station_sounding_obj.get_meanwind())
    print("Surface CAPE", station_sounding_obj.get_sb_cape())
    plot_stuve(station_sounding_obj, f"CM1 {name}")
    plot_skewT(station_sounding_obj, f"CM1 {name}")

if __name__ == "__main__":
    main()
