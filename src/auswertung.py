#!/usr/bin/python3

import matplotlib
matplotlib.use('Agg')
# Force matplotlib to not use any Xwindows backend.
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter,ScalarFormatter)
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from io import StringIO

import numpy as np
import numpy.ma as ma

from skewT import *
from readhtml_modul import *

import sys
#import datetime
import time
from copy import deepcopy


def datum ():
    hour = int(time.strftime("%H"))
    if hour >= 18:
        return time.strftime("%d.%m.%Y Time: 18 UTC")
    elif hour >= 12:
        return time.strftime("%d.%m.%Y Time: 12 UTC")
    elif hour >= 6:
        return time.strftime("%d.%m.%Y Time: 06 UTC")
    else:
        return time.strftime("%d.%m.%Y Time: 00 UTC")

#date=datum()

urlstring = "http://weather.uwyo.edu/cgi-bin/bufrraob.py?datetime=" + string1 + "%20" + string2 + ":00:00&id=" + string3 + "&type=TEXT:LIST"
