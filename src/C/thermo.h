/* 
This header include function to calculate thermodynamic parameter
of a sounding.
All function call back to a structure, 
which has before became the values of temperature, dewpoint, high for every press high.
The must function has as source sharpy.
*/

#ifndef THERMO_H
#define THERMO_H

/* *** Thermodynamische Funktionen *** */

/* Potential temperature of pressure level */

float theta (float *temp, float *press){ // einzelne Werte aus einer Struktur mit Call by Refrence methode.
/*  
    Parameters
    ----------
    temp : Temperature of a level or parcel (C)
    press : Pressure of a level or parcel (hPa)
    referrenz press : 1000 hPa
    Returns:
    --------
    pot : Potential temperature (K)
*/
    float pot;
    if(*press == 1000){
        return (*temp + K);
    } else {
        pot = (*temp + K) * pow((1000/(float) *press),ROCP);
        return pot;
    }
}

float pottemp (float temp, float druck){ // falls die Funktion mit call by value aufgerufen werden soll.
    
    return ( (temp+K) * pow((1000/druck),ROCP) );

}
/*
float theta (struct sounding_data *ptr){
    double temp;
    temp = ((*ptr).temp + K) * pow((1000/(double)(*ptr).druck),ROCP);
    return temp - K;
}
*/

float potlvl(float theta, float temp){
/*  

    Parameters
    ----------
    theta : Potential temperature of the parcel (K)
    temp : Temperature of the parcel (C)

    Returns
    -------
    Pressure Level (hPa [float]) of the parcel

    Funktionsweise
    --------------
    Die Formel für die Potentielle Temperatur 
    kann auch nach dem Druck aufgelöst werden:
    Cp (ln T - ln To) = R (ln p - ln po) 
    Cp/R ln(To/T) = po/p |To = Theata
    p= po / (To/T)^Cp/R
*/
    float thalvl;
    temp = temp + K;
/*
    Da bei trockenadiabatischer Erwärmung die Potentielle Temperatur erhalten ist, kann somit
    Theta als gegeben vorrausgesetzt werden und somit der Druck des lcl berechnet werden.
*/
    thalvl = 1000 / (pow((theta / temp),(1/ROCP))); 
    return thalvl;  
}

float drylift (float pottemp, float press){
/*
    Returns the temperature (C) of a parcel when raised to its LCL.

    Parameters
    ----------
    pottemp : Potentielle Temperatur (K)
    press : Pressure to which parcel is raised (hPa)

    Returns
    -------
    Temperature (K) of the parcel of drylift parcel at his new level.
*/
    float p0 = 1000.0;
    return (pottemp / pow((p0/press),ROCP) ); //-K
}

float lcltemp (float temp, float td){
/*
    Returns the temperature (C) of a parcel when raised to its LCL.

    Parameters
    ----------
    t : Temperature of the parcel (C)
    td : Dewpoint temperature of the parcel (C)

    Returns
    -------
    Temperature (C) of the parcel at it's LCL.
*/
    float dlt, s;// s = spread 
    s = temp - td;
    dlt = s * (1.2185 + 0.001278 * temp + s * (-0.00219 + 1.173e-5 * s -0.0000052 * temp)); //Quelle SHARPpy
    return temp - dlt;
}

float lclhigh (float t, float td, float high) {
    //return ((t-td) * 125); 
    return high + ((t-td) * 125.0);
}


/*void drylift(int druck, float temp, float td){
void drylift(struct sounding_data *struct_ptr){*/
/*  Lifts a parcel to the LCL and returns his new level and temperature.
    Parameters
    ----------
    pruck : Pressure of initial parcel in hPa
    temp : Temperature of inital parcel in C
    td : Dew Point of initial parcel in C
*//*
    Returns
    -------
    p2 : LCL pressure in hPa
    t2 : LCL Temperature in 

    (*struct_ptr).lclt = lcltemp((*struct_ptr).temp, (*struct_ptr).dewpoint);
    (*struct_ptr).lclp = potlvl(pottemp((*struct_ptr).temp, (float) (*struct_ptr).druck), (*struct_ptr).lclt);
    //printf("LCL(Lift Condensation Level):\tTemperatur:%.2f\tDruck:%.2f\n", ergebnis[0], ergebnis[1]);
    return;
}
*/

double wobf(double t){
/*
    Implementation of the Wobus Function for computing the moist adiabats.
    Parameters
    ----------
    t : Temperature (C)

    Returns
    -------
    Correction to theta (C) for calculation of saturated potential temperature.
*/  double npol, ppol;
    t = t - 20;
    if (t <= 0){
        npol = 1 + t * (-8.8416605e-3 + t * ( 1.4714143e-4 + t * (-9.671989000000001e-7 + t * (-3.2607217e-8 + t * (-3.8598073e-10)))));
        npol = 15.13 / (pow(npol,4));
        return npol;
    } else if (t > 0){
        ppol = t * (4.9618922e-07 + t * (-6.1059365e-09 + t * (3.9401551e-11 + t * (-1.2588129e-13 + t * (1.6688280e-16)))));
        ppol = 1 + t * (3.6182989e-03 + t * (-1.3603273e-05 + ppol));
        ppol = (29.93 / pow(ppol,4)) + (0.96 * t) - 14.8;
        return ppol;
    } else {
        return NAN;
    }
}

double satlift(double p, double thetam){
/*    
    Returns the temperature (C) of a saturated parcel (thm) when lifted to a
    new pressure level (hPa)

    Parameters
    ----------
    p : Pressure to which parcel is raised (hPa)
    thetam : Saturated Potential Temperature of parcel (C)

    Returns
    -------
    Temperature (C) of saturated parcel at new level
*/
    
    //if np.fabs(p - 1000.) - 0.001 <= 0: return thetam
    double eor = 999;
    double pwrp, rate, t1, t2, e1, e2;
    while (eor - 0.1 > 0){
        if (eor == 999) {                // First Pass
            pwrp = pow((p / 1000),ROCP);
            t1 = (thetam + 273.15) * pwrp - 273.15;
            e1 = wobf(t1) - wobf(thetam);
            rate = 1;
        }else{                         // Successive Passes
            rate = (t2 - t1) / (e2 - e1);
            t1 = t2;
            e1 = e2;
        }
        t2 = t1 - (e1 * rate);
        e2 = (t2 + 273.15) / pwrp - 273.15;
        e2 += wobf(t2) - wobf(e2) - thetam;
        eor = e2 * rate;
        if (eor < 0){
            eor *= -1;
        }
    }
    return t2 - eor;
}

double wetlift(float theta , float temp, float p2){
/*
    Lifts a parcel moist adiabatically to its new level.

    Parameters
    -----------
    press : Pressure of initial parcel (hPa). IF we have calculate the lclp we use the theta of the parcel, what until LCL is equal.
    So we can stint again the caculate of pottemp. theta : Pottemp bei given LCL.
    theta : Potential Temperature (K)
    temp : Temperature of initial parcel (C)
    p2 : Pressure of final level (hPa)
    Returns
    -------
    Temperature (C) by the pressure of final(new) level
*/
    double thetam;
    //theta = pottemp(t, p);
    thetam = (double) (theta -K) - wobf((double)(theta -K)) + wobf((double)temp);
    return (satlift((double) p2, (double) thetam));
   
}

double vappres(double temp){
    /*
    Returns the vapor pressure of dry air at given temperature

    Parameters  t : Temperature of the parcel (C)
    Returns ->  Vapor Pressure of dry air
    */
    double pol;
    pol = temp * (1.1112018e-17 + (temp * -3.0994571e-20));
    pol = temp * (2.1874425e-13 + (temp * (-1.789232e-15 + pol)));
    pol = temp * (4.3884180e-09 + (temp * (-2.988388e-11 + pol)));
    pol = temp * (7.8736169e-05 + (temp * (-6.111796e-07 + pol)));
    pol = 0.99999683 + (temp * (-9.082695e-03 + pol));
    return 6.1078 / pow(pol, 8);
}

double mixratio(float *press, float *temp){
    /*
    Returns the mixing ratio (g/kg) of a parcel

    Parameters
    ----------
    p : Pressure of the parcel (hPa)
    t : Temperature of the parcel (hPa)

    Returns
    -------
    Mixing Ratio (g/kg) of the given parcel

    */
    double x, wfw, fwesw, value;
    x = 0.02 * ((double) *temp - 12.5 + (7500. / (double) *press));
    wfw = 1. + (0.0000045 * (double) *press) + (0.0014 * x * x);
    fwesw = wfw * vappres((double) *temp);
    value = 621.97 * (fwesw / ((double) *press - fwesw));
    //printf("Sättigungsmischungsverhältnis in g/kg beträgt: %.2f\n", value);
    return value;
}

double virtemp(float *press, float *temp, float *td){
/*   
    Returns the virtual temperature (C) of a parcel.

    Parameters
    ----------
    p : The pressure of the parcel (hPa)
    t : Temperature of the parcel (C)
    td : Dew point of parcel (C)

    Returns:    Virtual temperature (C)
*/  double tk, w, vt;
    tk = *temp + K;
    w = 0.001 * mixratio(press, td);
    vt = (tk * (1 + w / EPS) / (1 + w)) - K;
    return vt;
}


float temp_at_mixrat( float w, float p){
    /*
    Returns the temperature (C) of air at the given mixing ratio (g/kg) and
    pressure (hPa)
    Parameters
    ----------
    w : number, numpy array
        Mixing Ratio (g/kg)
    p : number, numpy array
        Pressure (hPa)
    Returns
    -------
    Temperature (C) of air at given mixing ratio and pressure
    */
    double c1 = 0.0498646455;
    float c2 = 2.4082965 , c3 = 7.07475, c4 = 38.9114, c5 = 0.0915, c6 = 1.2035;
    float x = log10(w * p / (622.0 + w));
    x = (pow(10.0,((c1 * x) + c2)) - c3 + (c4 * pow((pow(10,(c5 * x)) - c6),2))) - K;
    return x;
}
/*Wet Bulb Temperature*/
float wet_bulb_temperature(float lclp, float lclt, float *press){
    /*
    Return the wet bulb temperature (C).
    
    Parameters
    ----------
    lclp : LCL Pressure for a given layer
    lclt : LCL Temperature for a given layer/level (C)
    press: Pressure of the initianal level
    Returns
    -------
    Wet Bulb Temperature (C)
    */
    float thta = NAN;
    thta = theta(&lclp, &lclt);
    return (float) wetlift(thta, lclt, *press);
}

/*Einige Thermodynamische Indizes*/

float lapserateneu (float high1, float temp1, float high2, float temp2){
    int highdiff=0;
    float tempdiff=NAN, lapse=NAN;
    
    if(high1 == high2){
        printf("False Input! \t Cannot calculate lapserate for one hight. \n");
        return NAN;
    }else if (high1 < high2){
        highdiff = high1 - high2;
        tempdiff = temp1 - temp2;
    }else if (high1 > high2){
        highdiff = high1 - high2;
        tempdiff = temp2 - temp1;
    }
    lapse = (tempdiff / (float) highdiff)*1000;
    return lapse;
}

void lapserate (struct sounding_data *ptr1, struct sounding_data *ptr2){ //Strukturen werden als Zeiger übegeben. Zugriff mit dem * Operator.
    int highdiff=0;
    float tempdiff=NAN;
    
    if(ptr1->high == (*ptr2).high){
        printf("False Input! %.1f and %.1f\n", ptr1->high, (*ptr2).high);
        (*ptr2).lapse = 0.0;
    }else if ((*ptr1).high < (*ptr2).high){ //(*ptr2).high übergibt Wert, prt2->high würde das gleich machen
        highdiff = (*ptr2).high - (*ptr1).high;
        tempdiff = (*ptr1).temp - (*ptr2).temp;
        (*ptr2).lapse = (tempdiff / (float) highdiff)*1000;
    }else if ((*ptr1).high > (*ptr2).high){
        highdiff = (*ptr1).high - (*ptr2).high;
        tempdiff = (*ptr2).temp - (*ptr1).temp;
        (*ptr2).lapse = (tempdiff / (float) highdiff)*1000;
    }
    //printf("%f", (*ptr2).lapse);
    return;
}


/* Brunt-Väisälä-Frequenz*/
/*
void brunt_frequnz (struct sounding_data *ptr1){
    
    (*ptr1).nz = sqrt( (G/(*ptr1.temp)) * ((*ptr2).lapse/10) -  0.96)

    return;
    }
*/

float k_index (float t_850, float t_500, float td_850, float t_700, float td_700){
    //Quelle: http://glossary.ametsoc.org
    return ( (t_850 - t_500) + td_850 - (t_700 - td_700) );
}

float koindex (float theat925, float theat850, float theat700, float theat500){
    //berechnung des ko-index: http://www2.wetter3.de/gfs025.html
    return (0.5*(theat700 + theat500 - (theat850 + theat925)) );
    //return 0.5*(theat700 + theat500 - (theat850 + theat1000) );
}

float tt (float t_850, float t_500, float td_850){
    return (t_850 - t_500) + (td_850 - t_500);
}

//hauptdruckflächen finden und funktionen übergeben
void research (struct sounding_data all[], struct indices *struct_ptr, int anzahl){
    int i=10;
    float theta925=NAN, high925=NAN, temp925=NAN;
    float theta850=NAN, high850=NAN, temp850=NAN, dew850=NAN;
    float theta700=NAN, high700=NAN, temp700=NAN, dew700=NAN;
    float theta500=NAN, high500=NAN, temp500=NAN;
    while (i < anzahl){
        if((int) all[i].druck == 925){
            theta925=all[i].thetae;
            high925=all[i].high;
            temp925=all[i].temp;
        }
        if((int) all[i].druck == 850){
            theta850=all[i].thetae;
            high850=all[i].high;
            temp850=all[i].temp;
            dew850=all[i].dewpoint;
        }   
        if((int) all[i].druck == 700){
            theta700=all[i].thetae;
            high700=all[i].high;
            temp700=all[i].temp;
            dew700=all[i].dewpoint;
        }
        if((int) all[i].druck == 500){
            theta500=all[i].thetae;
            high500=all[i].high;
            temp500=all[i].temp;
            break;
        }
        i++;
    }
    
    (*struct_ptr).lapse1 = lapserateneu(high500,temp500,high850,temp850);
    (*struct_ptr).lapse2 = lapserateneu(high500,temp500,high700,temp700);
    (*struct_ptr).lapse3 = lapserateneu(high850,temp850,high925,temp925);
    //(*struct_ptr).lapse4 = lapserateneu();
    //KO-Index
    (*struct_ptr).koindex = koindex (theta925, theta850, theta700, theta500);
    //K-Index
    (*struct_ptr).kindex = k_index(temp850,temp500,dew850,temp700,dew700);
    //The Total Totals index is attributable to Miller (1972)
    (*struct_ptr).tt = tt(temp850,temp500,dew850);
    return;
}

/*Integrationmittels Simsonsregel*/

double si_inetragtion (double a, double b){
    return (a + ((a+b)/(double) 2.0) + b)/(double)3.0;
}



/*CAPE*/


/*Luftpaket aufsteigen lassen*/
/*Wrapper für Pseudopotentielle Temperatur un CAPE*/


struct parcel * CAPE (struct parcel *cape,struct sounding_data drk[], float potvir, int index){
    /*
    This function calculate the convective avaialable potential energy
    Parameters:
    The function need the full sounding to check the energy for each level by each level
    Returns: 
    Save for each sounding level the Cape and Cin in the structure all[] and returns the structure mu-cape.
    mu-cape include sounding_datarmations of the must unstable cape (the level with the greatest energy or highest Theta-E values) 
    for example cin, eg, lcl, start level. 

    Parameter:
    start level : t(bzw. tv),td,press
    lcl : lclp
    struct sounding_data drk[anzahl];//Globale Variable
    index;
    modus : cape or capev
    */

    //struct for theoretical parcel
    
    int i=0, l=0, eq=0; //control variables (l for lfc and eq for equilibrium level)
    bool cin_flag = true; //Cin Flag
    double summecape = 0, summecin = 0; //variables to summation of each layer
    double b_bot=0, b_top=0, b=0; /*speed up*/
    float tv_parcel=0, td_parcel=0;/*temperature of the parcel*/
    //float tv_env=0; /*virtual temperature of the enviroment*/
    //float dp=0; /*Druck differenz*/

    i = index;
    b_bot = G * ( ( (double) ((*cape).temp)+K - ( (double) drk[i].virtemp+K) )/( (double) drk[i].virtemp+K)  );//Formel aus der Vorlesung Dynamik 1;
    i++;
    while( drk[i].druck > (*cape).lclp){ /*if pressure higher as lcl pressure we are still under the lcl*/
        
        tv_parcel = drylift(potvir,drk[i].druck);/*Trockenadiabatisch bis zum Druckhlevel in den Daten*/

        b_top = G * ( ( (double) (tv_parcel) - ( (double) drk[i].virtemp+K) )/( (double) drk[i].virtemp+K)  );
        b = ((b_bot + b_top)/(double) 2.0); /*Integration durch Trapezregel*/
        //printf("Speedup :%.2f\tLevel:%.1f\t%d\n",b,drk[i].druck,cin_flag);/*Output to controll*/
        if (b < 0){
            
            summecin += ( b * (double) (drk[i].high-drk[i-1].high)  );/*Zu Cin hinzu zählen*/
            //printf("%.1f\n",summecin);
            //cin_flag = true;
        }
        b_bot=b_top;
        i++;
    }

    /*Enviroment of the lcl*/
    /*wmean*/
    /*unkfioniert nicht,leider*/
/*

    dp = (float) (drk[i].druck - drk[i+1].druck);

    tv_env = ( (drk[i+1].virtemp+K)*((drk[i].druck-(*cape).lclp)/dp) ) + ( (drk[i].virtemp+K)*(((*cape).lclp - drk[i+1].druck)/dp) ); //K
    //printf("%.1f\t%.1f\n",((drk[i].druck-(*cape).lclp)/dp),(((*cape).lclp - drk[i+1].druck)/dp));
    //virtual temperature on parcel level
    tv_parcel = virtemp(&((*cape).lclp),&((*cape).lclt),&((*cape).lclt));
    //speedup im lcl
    b_top = G * ( ( (double) (tv_parcel + K) - (double) tv_env )/((double) tv_env)  );
    //printf("%d\t%.1f\t%.1f\t%.1f\n",i,tv_parcel,tv_env-K,b_top);
    //cin_flag seting
    if (b_top <= 0){ //if the speedup negativ in lcl set cin flag = 1
        cin_flag = true;
    }
    //Speedup for the area between the lcl and the datalevel under lcl
    b = ((b_bot + b_top)/(double) 2.0); //Integration durch Trapezregel

    if (b < 0){
        summecin += ( b * (double) ((*cape).lclh - drk[i].agl)  );/*Zu Cin hinzu zählen*/
        //lclh ist in agl gegeben. Für die Richtige Differenz muss hier entsprechend auch agl bei den daten verwendet werden.
        //cin_flag = true;
/*    }
    
    b_bot=b_top;

    //Speedup between next datalevel and lcl
    //Lift the parcel moist adiabatically of the parcel to the current level
    i++;

    td_parcel = (float) wetlift((*cape).pottemp, (*cape).lclt,drk[i].druck);
    tv_parcel = virtemp(&drk[i].druck,&td_parcel,&td_parcel);
    if (betrag( (tv_parcel+K) - (drk[i].virtemp+K)) > 10){
         warning("unrealistic values in CAPE");
    }
    
    b_top = G * ( ( (double) (tv_parcel + K) - ( (double) drk[i].virtemp+K) )/( (double) drk[i].virtemp+K)  );
    b = ((b_bot + b_top)/(double) 2.0); //Integration durch Trapezregel
    
    if ( (b <= 0) && (cin_flag == true) ){
        
        if (b < -3)
            warning("speedup is lower then -3 m/s^2");
        summecin += ( b * (double) (drk[i].agl - (*cape).lclh) ); //wenn negativ, dann die Beschleunigung zu CIN dazu rechnen
    } else if (b > 0) {
        /*if (cin_flag == true){
            cin_flag = false;
            l = i;
        }*/
/*        summecape += (b * (double) (drk[i].agl - (*cape).lclh) ); //Only the positiv values addition to the CAPE
        
        eq = i;
        if (b > 3)
            warning("speedup is greater then 3 m/s^2");

    } else if ( (b == NAN ) || (b > 10) ){
        printf("Fehler in CAPE!\n At level %f.1\n",drk[i].druck);
        summecape = NAN;
        summecin = NAN;
    }
    
    b_bot=b_top;


    i++;
*/
    /*Now, we calculate the parcel moist adiabatically to the next level until we catch 100 hPa*/
    while( (drk[i].druck >= 100) ){/*Estimate CAPE until 100 hPa.*/

        /*Lift the parcel moist adiabatically of the parcel to the current level*/
        td_parcel = (float) wetlift((*cape).pottemp, (*cape).lclt, drk[i].druck);
        tv_parcel = virtemp(&drk[i].druck,&td_parcel,&td_parcel);
        if (betrag( (tv_parcel+K) - (drk[i].virtemp+K)) > 30){
            //warning("unrealistic values in CAPE");
            //printf("%d\t%.1f\t%.1f\n",i,tv_parcel,drk[i].virtemp);
            break;
        }
        b_top = G * ( ( (double) (tv_parcel + K) - ( (double) drk[i].virtemp+K) )/( (double) drk[i].virtemp+K)  );/*Integrand*/

        b = ((b_bot + b_top)/(double) 2.0); /*Integration durch Trapezregel*/
        
        if ( (b <= 0) && (cin_flag == true) ){
            summecin += ( b * (double) (drk[i].high-drk[i-1].high)  ); //wenn negativ, dann die Beschleunigung zu CIN dazu rechnen
            if (b < -3)
                warning("speedup is lower then 3 m/s^2");

        } else if (b > 0) {
            if (cin_flag == true){
                cin_flag = false;
                l = i;
            }
            
            eq = i;
            summecape += (b * (double) (drk[i].high-drk[i-1].high) ); /*Only the positiv values addition to the */

            if (b > 3)
                warning("speedup is greater then 3 m/s^2");

        } else if ( (b == NAN ) || (b > 10)){
            printf("Fehler in CAPE!\n At level %f.1\n",drk[index].druck);
            summecape = NAN;
            summecin = NAN;
            break;/*Break because it is absurd*/
        }

        //printf("Speedup :%.2f\tLevel:%.1f\t%d\n",b,drk[i].druck,cin_flag);/*Output to controll*/
        b_bot=b_top;
        i++;//Set control variable to next layer
        
    }

    /*save the results*/
    if( summecape == 0) {
        (*cape).cape = 0;
        (*cape).cin = 0;
        (*cape).lfcp = 0; 
        (*cape).lfch = 0;
        (*cape).lfct = 0;
        (*cape).eqp = 0;
        (*cape).eqh = 0;
        (*cape).eqt = 0; 
    } else if ( summecape > 0){
        (*cape).cape = summecape;
        (*cape).cin = summecin;
        (*cape).lfcp = drk[l].druck; 
        (*cape).lfch = drk[l].agl;
        (*cape).lfct = drk[l].temp;
        (*cape).eqp = drk[eq].druck;
        (*cape).eqh = drk[eq].agl;
        (*cape).eqt = drk[eq].temp;   
    } else if ( summecape < 0 || summecape == NAN){
        printf("Fehler in CAPE!");
        (*cape).cape = NAN;
        (*cape).cin = NAN;
        (*cape).lfcp = NAN; 
        (*cape).lfch = NAN;
        (*cape).lfct = NAN;
        (*cape).eqp = NAN;
        (*cape).eqh = NAN;
        (*cape).eqt = NAN; 
    } 

    return cape;

}



struct parcel *mixed_layer(struct sounding_data all[],struct indices *severe){
    double potmean = 0;
    float mixmean = 0, druck=all[0].druck-50;
    short i=0;
    /*mean of Potential Temprature and mix ration*/
    while (all[i].druck > (all[0].druck - 100)){
        potmean += all[i].pottemp;
        mixmean += all[i].mix;
        i++;
    }
    
    (*severe).mmix = mixmean/ (float) i;
    //printf("%f\n%.1f\n",drylift((float) potmean/i,all[0].druck-50) - K,temp_at_mixrat(mixmean/i,all[0].druck-50));
    //printf("%f\t%.1f\n",potmean/i, tempmean);

    struct parcel *ml;
    ml = malloc(sizeof(struct parcel));
    
    (*ml).pottemp = potmean/ (float) i; //Arithmetisches Mittel der Potentiellen Temperatur
    (*ml).temp = drylift((*ml).pottemp,druck) - K; // Temperatur in der hälfte des mixed layers
    (*ml).dewp = temp_at_mixrat( mixmean/ (float) i,druck); //Taupunkt aus Mittel des mixing ratio
    (*ml).lclt = lcltemp((*ml).temp, (*ml).dewp); //LCL Temperatur berechnen
    (*ml).lclp = potlvl((*ml).pottemp, (*ml).lclt);
    (*ml).lclh = lclhigh((*ml).temp, (*ml).dewp,0);
    (*ml).thetae = (pottemp(wetlift((*ml).pottemp, (*ml).lclt,100),100));

    float potvir = virtemp(&druck, &(*ml).temp, &(*ml).dewp);//Kkurz noch Pottentielle Virtuelle Temperatur wird für Cin benötigt.
    ml = CAPE (ml,all,theta(&potvir,&druck),0);
    
    return ml;
}

int parcel_wrapper (struct sounding_data all[],struct parcel *mu, struct parcel *sb){
    int i=0;
    int ret = 0;
    struct parcel *run;
    //sb = malloc(sizeof(struct parcel));
    //mu = malloc(sizeof(struct parcel));
    run = malloc(sizeof(struct parcel));
    while (all[i].high < 3000){
    
        (*run).pottemp = all[i].pottemp;
        (*run).temp = all[i].temp;   
        (*run).dewp = all[i].dewpoint;    
        (*run).lclt = lcltemp(all[i].temp, all[i].dewpoint);
        (*run).lclp = potlvl(all[i].pottemp, (*run).lclt);
        (*run).lclh = lclhigh(all[i].temp, all[i].dewpoint,all[i].agl);
        (*run).thetae = (pottemp(wetlift((*run).pottemp, (*run).lclt,100),100));

        run = CAPE (run,all,all[i].potvir,i);

        all[i].thetae = (*run).thetae;
        all[i].cape = (*run).cape;
        all[i].cin = (*run).cin;
        //printf("%.1f\t",all[i].cin);
        if (i==0){
            (*sb).cape = (*run).cape;
            (*sb).cin = (*run).cin;
            (*sb).pottemp=(*run).pottemp;
            (*sb).thetae=(*run).thetae;
            (*sb).temp=(*run).temp;
            (*sb).dewp=(*run).dewp;
            (*sb).lclp=(*run).lclp;
            (*sb).lclh=(*run).lclh;
            (*sb).lclt=(*run).lclt;
            (*sb).lfcp=(*run).lfcp;
            (*sb).lfch=(*run).lfch;
            (*sb).lfct=(*run).lfct;
            (*sb).eqp=(*run).eqp;
            (*sb).eqh=(*run).eqh;
            (*sb).eqt=(*run).eqt;
        }
        if ( ((*mu).cape < (*run).cape)||(i==0) ){
            (*mu).cape = (*run).cape;
            (*mu).cin = (*run).cin;
            (*mu).pottemp=(*run).pottemp;
            (*mu).thetae=(*run).thetae;
            (*mu).temp=(*run).temp;
            (*mu).dewp=(*run).dewp;
            (*mu).lclp=(*run).lclp;
            (*mu).lclh=(*run).lclh;
            (*mu).lclt=(*run).lclt;
            (*mu).lfcp=(*run).lfcp;
            (*mu).lfch=(*run).lfch;
            (*mu).lfct=(*run).lfct;
            (*mu).eqp=(*run).eqp;
            (*mu).eqh=(*run).eqh;
            (*mu).eqt=(*run).eqt;
            ret = i;
        }
        i++;
    }    
    
    return ret;
}

float t_max (struct sounding_data all[]){
    float max_pot = 0;
    float superadiabatic = 2;
    short i=0;
    while (all[i].druck > (all[0].druck - 100)){
        if (max_pot < all[i].pottemp)
            max_pot = all[i].pottemp;
        i++;
    }
    return drylift(max_pot,all[0].druck) - K + superadiabatic;
}

#endif /*THERMO_H*/
