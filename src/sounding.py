#!/usr/bin/python3
import numpy as np
# project moduls
import src.meteolib as meteolib
from src.readhtml_bufr import readhtml2numpy
from src.utilitylib import load_yaml

# ---------------------------------------------------------------------------------------------------------------------


class sounding():
    def __init__(self, urlstring):
        config = load_yaml('config.yml', yaml_path='./src/')
        print(config)
        self.pres_temp, self.height, self.t_temp, self.dew_temp, self.mr_temp, winddir, self.wind_speed = readhtml2numpy(urlstring)
        self.agl = np.subtract(self.height, self.height[0])
        self.q_temp = meteolib.mixrat_to_q(self.mr_temp/1000)
        self.vir_temp = meteolib.virtuelle(self.q_temp, self.t_temp)
        self.thetae_temp = meteolib.thetae(self.pres_temp[:-1], self.t_temp, self.q_temp, p0=1000.)
        self.wetbulb_temp = self._calc_wetbulb()

        # kinematic transformations
        self.u_temp, self.v_temp = meteolib.uvwind(winddir, self.wind_speed)
        self.ushr = np.zeros(self.t_temp.size)
        self.ushr[1:] = np.subtract(self.u_temp[1:], self.u_temp[:-1])
        self.vshr = np.zeros(self.t_temp.size)
        self.vshr[1:] = np.subtract(self.v_temp[1:], self.v_temp[:-1])

        # strom relative

        # cape profile


    def _calc_surface_parcel(self):
        press_parcel = np.copy(self.pres_temp)
        lclp, lclt = meteolib.drylift(self.pres_temp[0], self.t_temp[0], self.dew_temp[0])
        mu2_parcel = np.zeros(press_parcel.size, dtype=np.float64)
        for i in range(0,press_parcel.size):
            if press_parcel[i] >= lclp:
                #mu_parcel[i] = thetas(theta(p_temp[int(sys.argv[2])], t_temp[int(sys.argv[2])], p2=1000.), press_parcel[i])
                mu2_parcel[i] = meteolib.thetas(meteolib.theta(lclp, lclt, p2=1000.), press_parcel[i])
            else:
                mu2_parcel[i] = meteolib.wetlift(lclp,lclt,press_parcel[i])


    # wetbulb
    def _calc_wetbulb(self):
        wetzero = np.zeros(np.size(self.t_temp))
        for i in range(0,np.size(self.t_temp)):
            wetzero[i] = meteolib.wetbulb(self.pres_temp[i], self.t_temp[i], self.dew_temp[i])
        return wetzero


    # mean wind
    def get_meanwind(self):
        return meteolib.mean_wind(self.u_temp, self.v_temp, self.pres_temp[:-1], stu=0, stv=0)


    # SRH
    def srh(self, sru, srv, high, mode=1):
        """
        Parameters:
        sru : storm relative wind u component
        srv : storm relative wind v component
        SRH = int_bot^top (ushr*srv - vshr*sru) dz
        """
        srh = 0
        h = 0
        i = 0

        while (self.agl[i] < high):
            if (mode == 1) or (mode == 2):
                srh += self.ushr[i] * ((srv[i] + srv[i+1])/2.0) - self.vshr[i] * ((sru[i] + sru[i+1])/2.0)
            elif (mode == 3):
                h = self.ushr[i] * ((srv[i] + srv[i+1])/2.0) - self.vshr[i] * ((sru[i] + sru[i+1])/2.0)
                if h >= 0:
                    srh += h
            else:
                srh += 0
            i += 1
            h = 0

        return srh


    def streamwise(self, sru, srv, high):
        h = 0
        i = 0
        H = 0
        while (self.agl[i] < high):
            h = self.ushr[i] * ((srv[i] + srv[i+1])/2.0) - self.vshr[i] * ((sru[i] + sru[i+1])/2.0)
            st = np.sqrt(((sru[i] +  sru[i+1])/2.0)**2 + ((srv[i] + srv[i+1])/2.0)**2)
            if (st != 0):
                H += h/st
            i += 1
            h = 0
        return H