#!/usr/bin/python3

import numpy as np

from io import StringIO
import sys

import urllib.request

def readhtml (urlstring,filestr,point_num=980):
    fp = urllib.request.urlopen(urlstring)
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()
    #sounding datas are between the PRE Tags
    #print(mystr)
    mystr = mystr.split("<PRE>")[1]
    mystr = mystr.split("</PRE>")[0][314:]

    #write the datas direktly in to the file
    fobj = open(filestr, 'w+')#open file
    fobj.write(mystr)
    fobj.close() #close file

    #pres, hght, temp, dwpt, relh, mixr, wdir, wspd, thta, thte, thtv = np.loadtxt(StringIO(mystr), delimiter=' ', unpack=True)

def readhtml2numpy (urlstring):
    """
    Parametrs:
    urlstring : str
    Return:
    tuple : (pressure,height,temperature,dewpoint,mixratio,wdir,wspd)
    This function convert the sounding data of wyoming to a numpy arrays and return them as a tuple
    """
    #read html file
    fp = urllib.request.urlopen(urlstring)
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()
    #sounding datas are between the PRE Tags
    #print(mystr)
    mystr = mystr.split("<PRE>")[1]
    mystr = mystr.split("</PRE>")[0][314:]
    #convert to a huge numpy array
    np.fromstring(mystr, dtype=float, sep='  ')
    #split numpy arra a part
    #pres 0 height 1 temp 2 dew 3 mix 5 wdir 6 wspd 7
    return (np.fromstring(mystr, dtype=float, sep='  ')[0::11],np.fromstring(mystr, dtype=float, sep='  ')[1::11],np.fromstring(mystr, dtype=float, sep='  ')[2::11],np.fromstring(mystr, dtype=float, sep='  ')[3::11],np.fromstring(mystr, dtype=float, sep='  ')[5::11],np.fromstring(mystr, dtype=float, sep='  ')[6::11],np.fromstring(mystr, dtype=float, sep='  ')[7::11])
"""
def get_data(point_num):
...     base_url = 'http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime=2020-04-04%200:00:00&id=10393&type=TEXT:LIST'
                    http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime=2020-04-04%200:00:00&id=10381&type=TEXT:LIST
...     #r = requests.get(base_url)
...     r = requests.get(base_url.format(point_num))
...     html_content = r.text
...     print(html_content)
... 
>>> get_data(4980)

"""
"""
import requests
def readhtml (urlstring,filestr,point_num=4980):
    r = requests.get(urlstring.format(point_num))
    mystr = r.text
    #mystr =  requests.get(urlstring.format(point_num)).text
    print(mystr)
    #sounding datas are between the PRE Tags
    #print(mystr)
    mystr = mystr.split("<PRE>")[1]
    mystr = mystr.split("</PRE>")[0][314:]

    #write the datas direktly in to the file
    fobj = open(filestr, 'w+')#open file
    fobj.write(mystr)
    fobj.close() #close file
"""
#http://weather.uwyo.edu/cgi-bin/bufrraob.py?src=bufr&datetime=2018-12-07%2000:00:00&id=10393&type=TEXT:LIST
#http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime=2020-04-04%200:00:00&id=10393&type=TEXT:LIST

urlstring = "http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime=" + sys.argv[1] + "%20" + sys.argv[2] + ":00:00&id=" + sys.argv[3] + "&type=TEXT:LIST"
#sys.argv[1]#2018-12-07
#sys.argv[2]#00
#sys.argv[3]#10393
#print(urlstring)
filestr = sys.argv[3] + ".txt"
readhtml(urlstring,filestr)



