#!/usr/bin/python3
import sys
from copy import deepcopy
import numpy as np

# project moduls
import src.meteolib as meteolib
from src.readhtml_bufr import readhtml2numpy
# ---------------------------------------------------------------------------------------------------------------------

def main():
    urlstring = (f"http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime={sys.argv[1]}%20{sys.argv[2]}:00:00"
                 f"&id={sys.argv[3]}&type=TEXT:LIST")
    pres, height, temp, dewpoint, mixing_ratio, winddir, wind_speed = readhtml2numpy(urlstring)
    spec_hum = mixrat_to_q(mixing_ratio/1000)
    vir_temp = virtuelle(spec_hum, temp)
    print(vir_temp)

if __name__ == "__main__":
    main()
