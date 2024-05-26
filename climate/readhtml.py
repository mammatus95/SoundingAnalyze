#!/usr/bin/python3

from io import StringIO
import urllib.request
import sys

#print(sys.argv)

def readhtml (urlstring,filestr,mode=1):
    fp = urllib.request.urlopen(urlstring)
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()

    seq = mystr.split("</PRE>")
    #print (seq[1])
    if ( mode==1):
        mystr = seq[0]
        #seq = mystr.split("")
        liste = seq[0].splitlines(True)
        fobj = open(filestr, 'w+')#öffnet Datei
        for i in range(10,len(liste)-2):
            if len(liste[i].split()) > 10:
                fobj.write(liste[i])
            #fobj.write("\n")
    
        fobj.close()#schließt Datei
    else:
        #print (seq[2])
        liste = seq[2].splitlines(True)
        fobj = open(filestr, 'w+')#öffnet Datei
        for i in range(6,len(liste)-1):
            if len(liste[i].split()) > 6:
                fobj.write(liste[i])
    
        fobj.close()#schließt Datei


#urlstring = "http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2018&MONTH=" + sys.argv[1] + "&FROM=" + sys.argv[2] +"00&TO=" + sys.argv[2] + sys.argv[3] + "&STNM=10548"

#print(urlstring)
# http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2020&MONTH=04&FROM=0412&TO=0412&STNM=10393
#http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2017&MONTH=08&FROM=2812&TO=2812&STNM=10238

if ( (sys.argv[4] != "10548") and (sys.argv[4] != "10410")):
    urlstring = "http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2018&MONTH=" + sys.argv[1] + "&FROM=" + sys.argv[2] + sys.argv[3] + "&TO=" + sys.argv[2] + sys.argv[3] + "&STNM=" + sys.argv[4]
    filestr = sys.argv[4] + ".txt"
    readhtml(urlstring,filestr)
else:
    filestr = sys.argv[4] + ".txt"
    ##braucht sonder behandlung
    if (int(sys.argv[3]) == 12):
        urlstring = "http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2018&MONTH=" + sys.argv[1] + "&FROM=" + sys.argv[2] +"00&TO=" + sys.argv[2] + sys.argv[3] + "&STNM=10548"
        readhtml(urlstring,filestr,2)
    else:
        urlstring = "http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2018&MONTH=" + sys.argv[1] + "&FROM=" + sys.argv[2] + sys.argv[3] + "&TO=" + sys.argv[2] + sys.argv[3] + "&STNM=" + sys.argv[4]
        readhtml(urlstring,filestr,1) 

    #print (urlstring)

