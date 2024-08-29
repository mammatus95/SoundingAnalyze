#!/usr/bin/python3
import numpy as np
import sys
import urllib.request

# ---------------------------------------------------------------------------------------------------------------------


def readhtml(urlstring, filestr):
    fp = urllib.request.urlopen(urlstring)
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()

    # sounding datas are between the PRE Tags
    mystr = mystr.split("<PRE>")[1]
    mystr = mystr.split("</PRE>")[0][314:]

    # write the datas direktly in to the file
    fobj = open(filestr, 'w+')  # open file
    fobj.write(mystr)
    fobj.close() # close file

# ---------------------------------------------------------------------------------------------------------------------


def readhtml2numpy(urlstring, dtype=np.float64):
    """
    Parametrs:
    urlstring : str
    Return:
    tuple : (pressure,height,temperature,dewpoint,mixratio,wdir,wspd)
    This function convert the sounding data of wyoming to a numpy arrays and return them as a tuple
    """
    # read html file
    fp = urllib.request.urlopen(urlstring)
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()
    # sounding datas are between the PRE Tags
    mystr = mystr.split("<PRE>")[1]
    mystr = mystr.split("</PRE>")[0][314:]
    mystr = mystr.replace("                                                              \n",
                        " -99.0  -99.0  -99.0  -99.0  -99.0  -99.0  -99.0  -99.0  -99.0\n")
    #print(mystr)
    # split numpy.array a part
    return (np.fromstring(mystr, dtype=dtype, sep='  ')[0::11],  # pressure
            np.fromstring(mystr, dtype=dtype, sep='  ')[1::11],  # height
            np.fromstring(mystr, dtype=dtype, sep='  ')[2::11],  # temperature
            np.fromstring(mystr, dtype=dtype, sep='  ')[3::11],  # dewpoint
            np.fromstring(mystr, dtype=dtype, sep='  ')[5::11],  # mixing ratio
            np.fromstring(mystr, dtype=dtype, sep='  ')[6::11],  # wind direction
            np.fromstring(mystr, dtype=dtype, sep='  ')[7::11])  # wind speed

# ---------------------------------------------------------------------------------------------------------------------
# sys.argv[1] : 2018-12-07
# sys.argv[2] : 00
# sys.argv[3] : 10393
# example links:
# http://weather.uwyo.edu/cgi-bin/bufrraob.py?src=bufr&datetime=2018-12-07%2000:00:00&id=10393&type=TEXT:LIST
# http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime=2020-04-04%200:00:00&id=10393&type=TEXT:LIST

def main():
    urlstring = (f"http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime={sys.argv[1]}%20{sys.argv[2]}:00:00"
                 f"&id={sys.argv[3]}&type=TEXT:LIST")
    print(urlstring)
    file_str = f"{sys.argv[3]}.txt"
    readhtml(urlstring, file_str)
    pres, height, temp, dewpoint, mixing_ratio, winddir, wind_speed = readhtml2numpy(urlstring)
    print(pres.max(), height.size, temp.size, dewpoint.size, mixing_ratio.size, winddir.size, wind_speed.size)


if __name__ == "__main__":
    main()
