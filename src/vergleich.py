#!/usr/bin/python3

import time
import sys
import urllib.request
import os
import matplotlib
matplotlib.use('Agg')
# Force matplotlib to not use any Xwindows backend.
from matplotlib import pyplot as plt
import numpy as np
import numpy.ma as ma

#class temp:#
#
#    def __init__ (self, t):
#        self.name = name
####Trockenadiabte Plotten####
ROCP = 0.2855721
K = 273.15

def potlvl ( theta, temp):
    """
    Returns the level (hPa) of a parcel.

    Parameters
    ----------
    theta : Potential temperature of the parcel (C)
    temp : Temperature of the parcel (C)

    Returns
    -------
    Pressure Level (hPa [float]) of the parcel
    """

    temp = temp + K;
    theta = theta + K;
    thalvl = 1000.0 / (pow((theta / temp),(1.0/ROCP))); 
    return thalvl;
#################################

def vergleich (urlstring,string):
    fp = urllib.request.urlopen(urlstring)
    mybytes = fp.readlines()
    mystr=[]
    for i in range(1,len(mybytes)-1):
        mystr.append(mybytes[i].decode("utf8"))
        #print(mybytes[i].decode("utf8"))
    fp.close()

    i = 0
    filestr=[]
    title=[]
    legende=[]
    dateiname=[]
    for k in range(0,3):
        j = mystr.index("<PRE>\n",i)
        n = mystr.index("</PRE>\n",i)

        name = mystr[j-1].split()
        #print(name)

        #verschiedene string die später noch gebraucht werden
        
        filestr.append(name[0][4:9] + "_" + name[len(name)-4][0:2] + ".txt")
        dateiname.append("cmp" + name[0][4:9] + name[len(name)-4][0:2] + name[len(name)-3] + name[len(name)-2] + name[len(name)-1][0:4] + ".png")
        title.append(name[0][4:9])
        legende.append("kurve vom " + name[len(name)-3] + "." + name[len(name)-2] + "  " + name[len(name)-4][0:2] + " UTC")

        fobj = open(filestr[k], 'w+')
    
        for i in range(j+5,n):
            if len(mystr[i].split()) > 8:
                fobj.write(mystr[i])
        fobj.close()
        i=n+1
    return filestr, dateiname, title, legende

#Current date + time
hour = int(time.strftime("%H"))
day = int(time.strftime("%d"))
month = int(time.strftime("%m"))
year = int(time.strftime("%Y"))

if (hour < 13 and day == 1):
    start = str(day-1) + "12"
    end = str(day) + "00"
    month = str(month-1)
    year = str(year)
elif (hour < 13 ):
    start = str(day-1) + "12"
    end = str(day) + "00"
    month = str(month)
    year = str(year)
elif (hour >= 13):
    start = str(day) + "00"
    end = str(day+1) + "00"
    month = str(month)
    year = str(year)
#fertig
#station="11120"#"10393"
station=sys.argv[1]
#Url String erzeugen:
#urlstring = "http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2017&MONTH=08&FROM=2712&TO=2800&STNM=10393"
urlstring = "http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR="+year+"&MONTH=" + month + "&FROM=" + start + "&TO=" + end + "&STNM=" + station
#print(urlstring)
#auseinander schnippel
filestr, dateiname, title, legende = vergleich(urlstring,station)
print(filestr, dateiname, title, legende)
#weiter gehts
#C ausführen
for k in range(0,len(filestr)-1):
    name = "./vergleich " + filestr[k] + " " + filestr[len(filestr)-1]
    os.system(name)

    #ergebnis laden
    high1, druck1, temp1, td1, high2, druck2, temp2, td2 = np.loadtxt("vergleich.txt",delimiter=',', unpack=True)
    #print(data)
    #Plotten
    start = legende[k]#gleich elmente wie bei C call
    end = legende[len(filestr)-1]

    b, ax = plt.subplots(1,figsize=(15,11))
    plt.grid(linewidth=1.9)

    ax.plot(temp1,druck1,'-g',td1,druck1,'--g', temp2,druck2,'-m',td2,druck2,'-.m',lw=2)
    ax.legend(( "Temperatur" + start, "Taupunkts" + start , "Temperatur" + end, "Taupunkts" + end),loc=1)
    ax.set_title(title[len(filestr)-1] + ": Vergleich der Aufstiege von 12 und 00 UTC")
    plt.yscale('log')
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks(np.linspace(100,1000,10))
    ax.set_ylim(1050,100)
    ax.set_xticks(np.arange(-90,40,10))
    ax.set_xlim(-90,40)
    plt.xlabel("Temperatur [C]")
    plt.ylabel("Druck [hPa]")

    #trockenadiabten plotten
    pot = np.arange(-80,150,10)
    for i in range(0,np.size(pot)):
        temp = np.arange(-90,pot[i]+20,10)
        pres = np.ones(np.size(temp))
        for j in range(0,np.size(temp)):
            pres[j] = potlvl(pot[i],temp[j])
        ax.plot(temp,pres,'-g',lw=0.5)


    #plt.show()
    #plt.savefig(dateiname[2])
    plt.savefig(title[2] + "_" + str(k) +"ver.png")
    plt.close()
