#!/usr/bin/python3
import numpy as np
# project moduls
import src.meteolib as meteolib
from src.readhtml_bufr import readhtml2numpy
from src.utilitylib import load_yaml

# ---------------------------------------------------------------------------------------------------------------------


class sounding():
    # env : sounding profile
    # parcel : 
    def __init__(self, urlstring):
        config = load_yaml('config.yml', yaml_path='./src/')
        print(config)
        self.pres_env, self.height, self.t_env, self.dew_env, self.mr_env, winddir, self.wind_speed = readhtml2numpy(urlstring)
        index_array = np.where(self.pres_env > config["pmin"])
        self.pres_env = self.pres_env[index_array]
        self.height = self.height[index_array]
        self.t_env = self.t_env[index_array]
        self.dew_env = self.dew_env[index_array]
        self.mr_env = self.mr_env[index_array]
        winddir = winddir[index_array]
        self.wind_speed = self.wind_speed[index_array]

        self.mr_env /= 1000
        self.t_env += meteolib.cr['ZEROCNK']
        self.dew_env += meteolib.cr['ZEROCNK']

        self.agl = np.subtract(self.height, self.height[0])
        self.q_env = meteolib.mixrat_to_q(self.mr_env)
        self.tvir_env = meteolib.virtuelle(self.q_env, self.t_env)
        self.theta_env = meteolib.theta(self.pres_env, self.t_env, p0=1000.)
        self.thetae_env = meteolib.thetae(self.pres_env, self.t_env, self.q_env, p0=1000.)
        self.wetbulb_env = self._calc_wetbulb()

        # kinematic transformations
        self.u_env, self.v_env = meteolib.uvwind(winddir, self.wind_speed)
        self.ushr = np.zeros(self.t_env.size)
        self.ushr[1:] = np.subtract(self.u_env[1:], self.u_env[:-1])
        self.vshr = np.zeros(self.t_env.size)
        self.vshr[1:] = np.subtract(self.v_env[1:], self.v_env[:-1])

        # strom relative

        # sb cape
        self.SB_CAPE, self.SB_CIN, self.SB_LFC, self.SB_EL = self._compute_lift_parcel(0, self.pres_env[0], self.t_env[0], self.dew_env[0])

        # cape profile


    def _compute_lift_parcel(self, index_start, start_pres, start_t, start_dew):
        press_parcel = np.copy(self.pres_env[index_start:])
        lclp, lclt = meteolib.drylift(start_pres, start_t, start_dew)
        print(lclp, lclt)
        tv_parcel = np.zeros(press_parcel.size, dtype=np.float64)
        for i in range(0, press_parcel.size):
            if press_parcel[i] >= lclp:  # unsaturated 
                tv_parcel[i] = meteolib.thetas(meteolib.theta(lclp, lclt, p0=1000.), press_parcel[i])
            else:  # saturated 
                tv_parcel[i] = meteolib.wetlift(lclp, lclt, press_parcel[i])  # not virtuel temperature

        buoyancy = self._compute_buoyancy(index_start, press_parcel, tv_parcel)
        print(buoyancy.max(),buoyancy.min())
        CAPE, CIN, LFC, EL = self._compute_CAPE(index_start, buoyancy)
        return CAPE, CIN, LFC, EL



    def _compute_buoyancy(self, index_start, press_parcel, tv_parcel):
        buoyancy = np.zeros(press_parcel.size, dtype=np.float64)
        #i = 0
        #while(press_parcel[i] >= 100):
        #    buoyancy[i] = meteolib.cr['G'] * ((tv_parcel[i] - self.tvir_env[i+index_start])/self.tvir_env[i+index_start])
        #    i += 1
        buoyancy = meteolib.cr['G'] * np.divide(np.subtract(tv_parcel, self.tvir_env[index_start:]), self.tvir_env[index_start:])
        return buoyancy

    def _compute_CAPE(self, index_start, buoyancy):
        dz = np.subtract(self.height[index_start+1:], self.height[index_start:-1])
        CAPE = np.nan
        CIN = np.nan
        LFC = np.nan
        EL = np.nan
        if np.nanmax(buoyancy) > 0:
            # CAPE will be the total integrated positive buoyancy
            B_pos = np.copy(buoyancy)
            B_pos[np.where(B_pos < 0)] = 0
            CAPE = np.nansum(np.multiply(dz, 0.5*np.add(B_pos[1:], B_pos[:-1])))
            # CIN will be the total negative buoyancy below the height of maximum
            B_max = np.nanmax(buoyancy)
            idx_max = np.where(buoyancy == B_max)[0][0]
            B_neg = np.copy(buoyancy)
            B_neg[np.where(B_neg > 0)] = 0
            B_neg[idx_max:] = 0
            CIN = np.nansum(np.multiply(dz, 0.5*np.add(B_pos[1:], B_pos[:-1])))

            # LFC will be the last instance of negative buoyancy before the
            # continuous interval that contains the maximum in buoyancy
            fneg = np.where(buoyancy < 0)[0]
            inn = np.where(fneg < idx_max)[0]
            fneg = fneg[inn]

            if len(fneg) > 0:
                LFC = 0.5*(self.height[np.max(fneg)] + self.height[np.max(fneg)+1])
            else:
                LFC = self.height[index_start]

            # EL will be last instance of positive buoyancy
            fpos = np.where(buoyancy > 0)[-1]
            EL = 0.5*(self.height[np.max(fpos)] + self.height[np.max(fpos)+1])
        else:
            CAPE = 0
            CIN = 0
            LFC = np.nan
            EL = np.nan

        return CAPE, CIN, LFC, EL

    # wetbulb
    def _calc_wetbulb(self):
        wetzero = np.zeros(np.size(self.t_env))
        for i in range(0,np.size(self.t_env)):
            wetzero[i] = meteolib.wetbulb(self.pres_env[i], self.t_env[i], self.dew_env[i])
        print(wetzero)
        return wetzero


    # mean wind
    def get_meanwind(self):
        return meteolib.mean_wind(self.u_env, self.v_env, self.pres_env, stu=0, stv=0)

    # surface based CAPE
    def get_sb_cape(self):
        return self.SB_CAPE

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