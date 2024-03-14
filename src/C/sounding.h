/*   sounding.h      */

#ifndef SOUNDING_H
#define SOUNDING_H


#include <math.h>
#include <stdbool.h>//für cin_flag
#include <time.h>

/*Konstanten*/
//#define ROCP 0.28571426
#define ROCP 0.2855721 //287/1005
//0.2858566 
#define K 273.15
#define EPS 0.62197
#define PI 3.14159265
#define G 9.80665 //m*s^(-2)

/*Makros*/
#define SQR(a) ((a)*(a))
#define KT2MS(v) (v*0.514444)


//structure with a values of a pressure level
struct sounding_data{
    float druck;    //hPa
    float high;     //m
    float agl;      //m
    float temp;     //C
    float dewpoint; //C
    double mix;     //g/kg
    float pottemp;  //K
    float thetae;   //K
    float virtemp;  //C
    float potvir;   //K
    float wetbulb;  //C
    float lapse;    //K/km
    float cape;     //J/kg 
    float cin;      //J/kg
    float wdir;     //m/s
    float wspd;     //Degres
    float u;        //m/s
    float v;        //m/s
    float ushr;     //m/s
    float vshr;     //m/s
    float srul;     //m/s
    float srvl;     //m/s
    float srur;     //m/s
    float srvr;     //m/s
};

struct stormmotion {
    float rstu; //m/s
    float rstv; //m/s
    float lstu; //m/s
    float lstv; //m/s
    float sr_shear1;  //m/s
    float sr_shear3;  //m/s
    float sr_shear6;  //m/s
    float sr_shear24; //m/s
    float srh1; //right //m2/s2
    float srh2; //right
    float srh3; //right
    float srh1p; //right
    float srh2p; //right
    float srh3p; //right
    float leftsrh1;
    float leftsrh2;
    float leftsrh3;
    float stw1;
    float stw2;
    float stw3;
    float angle;
};

struct parcel {
    float pottemp; //Potential temperature of the start level
    float thetae;  //pseudo-potentielle temperature of the parcel
    float temp;    //temperature of the start level
    float dewp;    //dewpoint of the start level
    float cape;    //convective available potential energy
    float cin;     //convective inhibition
    float lclp;    //lifting condensation level pressure
    float lclh;    //lifting condensation level high
    float lclt;    //lifting condensation level temperature
    float lfcp;    //level of free convection pressure
    float lfch;    //level of free convection high
    float lfct;    //level of free convection temperature
    float eqp;     //equilibrium level pressure
    float eqh;     //equilibrium level hight
    float eqt;     //equilibrium level temperature
};

struct inflow {
    float ptop; //hPa
    float pbot; //hPa
    float esrh; //m2/s2
    float ebwd; //m/s
};

struct indices {
    /*Allgmeines*/
    //float wmean_u; //pressure-weighted mean wind 0-6 km
    //
    /*Thermische Indizes*/
    struct parcel sb;
    struct parcel mu;
    struct parcel ml;
    float wetbulbzero; //m
    float koindex;
    float kindex;
    float tt;
    float lapse1; // K/km //lapserateneu(high500,temp500,high850,temp850);
    float lapse2; //lapserateneu(high500,temp500,high700,temp700);
    float lapse3; //lapserateneu(high850,temp850,high925,temp925);
    float tx;
    float mmix; //g/kg
    /*Shear*/
    struct stormmotion bunker;
    float hel500m; //m2/s2
    float hel1;
    float hel2;
    float hel3;
    float shear500m; //m/s
    float shear1;
    float shear3;
    float shear6;
    float shear8;
    float length1;
    float length28;
    float length3;
    float length6;
    float length8;
    /*Severe weather indices*/
    float ehi;
    float wxshr;
    float brn_shear; //m2/s2
    float brn_mu;
    float brn_ml;
    float ship;
    struct inflow efi;
    float scp;
    float stp;
};

/*Allgemeine Funktionen*/

//Umrechnung von Knoten in m/s; Das Makro dazu ist unter Makros zufinden
/*float kt2ms(float v){
    return v * 0.514444;
}*/

float betrag(float x){
    if (x < 0){
        return x * (-1);
    }else if (x >= 0){
        return x;
    } else {
        return NAN;
    }
}

void warning (char *string) {

    fprintf(stderr, "Warning : %s\n",string); //stderr Standardfehlerausgabe

    return;
}

/*local headers*/

#include "thermo.h"
#include "shear.h"
#include "ausgabe.h"

/*some more indicies*/

//temp.ehi = ehi(temp.sharppy.srh3p, temp.mu.cape);
float ehi ( float hel, float cape){
    return ( (hel * cape) / (float) 160000);
}

//temp.ship = ship (drk, temp.shear6, &temp.mu, temp.lapse2,j);
float ship ( struct sounding_data drk[], float shr6, int ret, float *mucape, float laps, int j){
    /*
    Calculate the Sig Hail Parameter (SHIP) by SPC

    Parameters
    ----------
    
    mucape : (optional) Most-Unstable Parcel
    dw : dewpoint of must unstable parcel at beginn of lifting
    laps : (optional) 700 - 500 mb lapse rate (C/km)
    t500 : (optional) 500 mb temperature (C)
    shr6 : (optional) 0-6 km shear (m/s)
    frz_lvl : (optional) freezing level (m)
    */
    int i=0;
    float t500=0,dw=11, frz_lvl = NAN, ship=0;
    while (i<j){
        if ( drk[i].druck == 500)
            t500 = drk[i].temp;
        if ( drk[i].temp <= 0 )
            frz_lvl = drk[i].high;
        i++;
    }
    dw = drk[ret].dewpoint;

    if (shr6 > 27.0)
        shr6 = 27.0;
    if (shr6 < 7.0)
        shr6 = 7.0;

    if (dw > 13.6)
        dw = 13.6;
    if (dw < 11.0)
        dw = 11.0;

    if (t500 > -5.5)
        t500 = -5.5;

    ship = (-1.0) * ( ( *mucape * dw * laps * t500 * shr6) / (float) 42000000 );
    
    if ( *mucape < 1300)
        ship = ship*(*mucape/ (float) 1300);
    
    if (laps < 5.8)
        ship = ship*( laps / (float) 5.8);

    if (frz_lvl < 2400.0)
        ship = ship * (frz_lvl/ (float) 2400);

    return ship;
}


//temp.scp = scp (&temp.mu.cape, temp.efi.esrh, temp.efi.ebwd);
float scp (float *mucape, float esrh, float ebwd){
    /*
    SCP = (muCAPE / 1000 J/kg) * (ESRH / 50 m2/s2) * (EBWD / 20 m/s)
    additionally:
    - EBWD is divided by 20 m s-1 in the range of 10-20 m s-1. 
    - EBWD less than 10 m s-1 is set to zero.
    - EBWD greater than 20 m s-1 is set to one.
    */

    if ( (ebwd < 10 ) || (ebwd == NAN ) )
        ebwd = 0;
    else if ((ebwd >= 10) && (ebwd <= 20))
        ebwd /= (float) 20;
    else if (ebwd > 20 )
        ebwd = 1;
    
    return (*mucape / (float) 1000) * ( esrh / (float) 50 ) * ebwd;
}
/*
float stp_cin (){

//Sig Tor (CIN) = (mlCAPE / 1500 J/kg) * (ESRH / 150 m2/s2) * (EBWD / 12 m/s) * ((2000 - mlLCL) / 1000) * ((mlCINH + 200) / 150)


}

float stp_fixed (){


}
*/

float wxshear (float *mlcape, float *shear6){
    /*
      ML WMAX * SHEAR
      ML WMAX = sqrt(2*CAPE) * Bulk Shear 0-6 km
    */
    if (*mlcape > 0){
        return (sqrt(2 * (*mlcape)) * (*shear6));
    } else {
        return 0;
    }
}

void brnshear (struct sounding_data drk[], struct indices *severe){
    /*
    Calculation of the bulk Richardson number (BRN), which is used to 
    forecast storm type. We obtain BRN by division of CAPE through the BRN Shear.
    The BRN Shear is the vector difference between
    the 0-6km mean wind and the 0-500m mean wind.   
    */

    float du=0,dv=0;
    // Calculate the lowest 500m mean wind
    du = mean_wind_u(drk,0,500)-mean_wind_u(drk,0,6000);
    dv = mean_wind_v(drk,0,500)-mean_wind_v(drk,0,6000);

    (*severe).brn_shear = SQR(sqrt( SQR(du) + SQR(dv) ) ) / 2.0;
    (*severe).brn_mu   = (*severe).mu.cape/(*severe).brn_shear;
    (*severe).brn_ml   = (*severe).ml.cape/(*severe).brn_shear;
    return;
}

//starten der einzelen Funktionen zum Auswerten der Daten
void calculate(struct sounding_data *struct_ptr){
    float lclt=0,lclp=0;
    uvwind(struct_ptr);
    (*struct_ptr).pottemp = theta(&(*struct_ptr).temp,&(*struct_ptr).druck);//call by refrence
    (*struct_ptr).virtemp = virtemp(&(*struct_ptr).druck, &(*struct_ptr).temp, &(*struct_ptr).dewpoint);//call by refrence
    (*struct_ptr).potvir = theta(&(*struct_ptr).virtemp,&(*struct_ptr).druck);//Potentielle Temperatur über die Virtuelle Temperatur
    (*struct_ptr).mix = mixratio(&(*struct_ptr).druck, &(*struct_ptr).dewpoint);//call by refrence
    lclt = lcltemp((*struct_ptr).temp, (*struct_ptr).dewpoint);
    lclp = potlvl((*struct_ptr).pottemp, lclt);
    /*(*struct_ptr).lclh = lclhigh((*struct_ptr).temp, (*struct_ptr).dewpoint, (*struct_ptr).high);*//*Höhe beim Must unstable Cape zu definieren macht kein Sinn*/
    (*struct_ptr).thetae = (pottemp(wetlift((*struct_ptr).pottemp, lcltemp((*struct_ptr).temp, (*struct_ptr).dewpoint),100),100));
    wet_bulb_temperature(lclp,lclt,&(*struct_ptr).druck);
    return;
}

#endif /*SOUNDING_H*/
