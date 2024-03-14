#ifndef AUSGABE_H
#define AUSGABE_H

void print_parcel(struct parcel *ptr){
    printf("convective available potential energy :\t%.1f\n",(*ptr).cape);
    printf("convective inhibition :\t%.1f\n",(*ptr).cin);
    printf("Potential temperature of the start level :\t%.1f\n",(*ptr).pottemp);
    printf("pseudo-potentielle temperature of the parcel :\t%.1f\n",(*ptr).thetae);
    printf("temperature of the start level :\t%.1f\n",(*ptr).temp);
    printf("dewpoint of the start level :\t%.1f\n",(*ptr).dewp);
    printf("lifting condensation level pressure :\t%.1f\n",(*ptr).lclp);
    printf("lifting condensation level high :\t%.1f\n",(*ptr).lclh);
    printf("lifting condensation level temperature :\t%.1f\n",(*ptr).lclt);
    printf("level of free convection pressure :\t%.1f\n",(*ptr).lfcp);
    printf("level of free convection high :\t%.1f\n",(*ptr).lfch);
    printf("level of free convection temperature :\t%.1f\n",(*ptr).lfct);
    printf("equilibrium level pressure :\t%.1f\n",(*ptr).eqp);
    printf("equilibrium level hight :\t%.1f\n",(*ptr).eqh);
    printf("equilibrium level temperature :\t%.1f\n",(*ptr).eqt);
    printf("\n");
}

//Ausgabe der berechneten Werte
void ausgabe1 (struct sounding_data *struct_ptr, int i){
    if(i == 0){
        printf("\n\tRichtung Geschwindigkeit\tu \t v\tushr\tvshr\tsru\tsrv\tcin\tcape\n");
        printf("\t\t\t m/s \t\t m/s \t m/s \tm/s \tm/s \tm/s \tm/s \t J/kg\t J/kg\n");
        printf("\t----------------------------------------------------------------------------------------------\n");
    }

/*    printf("\t%f", (*struct_ptr).druck);
    printf("\t%f", (*struct_ptr).high);
    printf("\t%.1f",(*struct_ptr).dewpoint);
    printf("\t%.1f", (*struct_ptr).temp);
    printf("\t%.1f",(*struct_ptr).virtemp);
    printf("\t%.1f",(*struct_ptr).lclp);
    printf("\t%.1f",(*struct_ptr).lclt);
    printf("\t  %.1f",(*struct_ptr).pottemp);
    printf("\t\t%.1f",(*struct_ptr).potvir);
    printf("\t  %.1f   ",(*struct_ptr).thetae);
    printf("%.1f\t",(*struct_ptr).lapse);*/
    printf("\t%.1f\t\t",(*struct_ptr).wdir);
    printf("%.1f\t\t",(*struct_ptr).wspd);
    printf("%.1f\t",(*struct_ptr).u);
    printf("%.1f ",(*struct_ptr).v);
    printf("\t%.1f",(*struct_ptr).ushr);
    printf("\t%.1f",(*struct_ptr).vshr);
    printf("\t%.1f",(*struct_ptr).srur);
    printf("\t%.1f",(*struct_ptr).srvr);
//    printf("\t%.1f",(*struct_ptr).srul);
//    printf("\t%.1f",(*struct_ptr).srvl);
    printf("\t%.1f",(*struct_ptr).cin); 
    printf("\t%.1f\n",(*struct_ptr).cape); 
    return;
}


//Ausgabe in die Datei radiodata.txt
void fileoutput (struct sounding_data *struct_ptr, FILE *data){

    fprintf(data,"%7.1f,", (*struct_ptr).agl);
    fprintf(data,"%7.2f,", (*struct_ptr).druck);
    fprintf(data,"%.1f,", (*struct_ptr).temp);
    fprintf(data,"%.1f,",(*struct_ptr).dewpoint);
    fprintf(data,"%.1f,",(float) (*struct_ptr).mix);
    fprintf(data,"%.1f,", (*struct_ptr).virtemp);
    fprintf(data,"%.1f,",((*struct_ptr).pottemp) - K);
    fprintf(data,"%.1f,",((*struct_ptr).thetae) - K);
    fprintf(data,"%.1f,",(*struct_ptr).lapse);
    fprintf(data,"%.1f,",(*struct_ptr).cape);
    fprintf(data,"%.1f,",(*struct_ptr).cin);
    fprintf(data,"%.1f,", (*struct_ptr).wspd);
    fprintf(data,"%.1f,", (*struct_ptr).u);
    fprintf(data,"%.1f,", (*struct_ptr).v);
    fprintf(data,"%.1f,", (*struct_ptr).srur);
    fprintf(data,"%.1f,", (*struct_ptr).srvr);
    fprintf(data,"%.1f,", (*struct_ptr).ushr);
    fprintf(data,"%.1f,", (*struct_ptr).vshr); 
    fprintf(data,"%.1f,", (*struct_ptr).srul);
    fprintf(data,"%.1f\n", (*struct_ptr).srvl);
    return;
}

//Ausgabe der Parameter/Indizes
void ausgabe2 (struct indices *ptr, struct tm *t){
    printf("Monat: %d Tag: %d Uhrzeit: %d UTC\n", ((*t).tm_mon + 1), (*t).tm_mday,(*t).tm_hour);
    printf("Thermodynamsiche Inizes\n");
    printf("Lapse Rate zwischen 850 & 500\t\t:%.1f\tK/KM\n", (*ptr).lapse1);
    printf("Lapse Rate zwischen 700 & 500\t\t:%.1f\tK/KM\n", (*ptr).lapse2);
    printf("Lapse Rate zwischen 925 & 850\t\t:%.1f\tK/KM\n", (*ptr).lapse3);    
    printf("KO-Index:\t:%.1f\tK\n", (*ptr).koindex);
    printf("SHIP\t:%.1f\t\n", (*ptr).ship);
    printf("T Max\t:%.1f\tC\n\n", (*ptr).tx);
    printf("Mean Mix Ratio\t:%.1f\tC\n\n", (*ptr).mmix);

    printf("                         Cape  CIN  J/kg  High of : LCL  LFC  EQ\n");
    printf("Surface base Parcel  : %.1f  %.1f\t\t%.1f\t%.1f\t%.1f\n",(*ptr).sb.cape,(*ptr).sb.cin,(*ptr).sb.lclh,(*ptr).sb.lfch,(*ptr).sb.eqh);
    printf("Most Unstable Parcel : %.1f  %.1f\t\t%.1f\t%.1f\t%.1f\n",(*ptr).mu.cape,(*ptr).mu.cin,(*ptr).mu.lclh,(*ptr).mu.lfch,(*ptr).mu.eqh);
    printf("Mixed Layer Parcel   : %.1f  %.1f\t\t%.1f\t%.1f\t%.1f\n\n",(*ptr).ml.cape,(*ptr).ml.cin,(*ptr).ml.lclh,(*ptr).ml.lfch,(*ptr).ml.eqh);

    printf("\nWind\n");
    printf("Bulk Shear:\n");
    printf("Shear 0-1km\t:%.1f\tm/s\n",(*ptr).shear1);
    printf("Shear 0-3km\t:%.1f\tm/s\n",(*ptr).shear3);
    printf("Shear 0-6km\t:%.1f\tm/s\n",(*ptr).shear6);
    printf("Helicity0-1km\t:%.1f\tm2/s2\n", (*ptr).hel1);
    printf("Helicity0-2km\t:%.1f\tm2/s2\n", (*ptr).hel2);
    printf("Helicity0-3km\t:%.1f\tm2/s2\n", (*ptr).hel3);
    printf("Energy-Helicity Index\t:%.1f\tm2/s2\n", (*ptr).ehi);
    printf("Effektive Inflow Layer\n");
    printf("\tESRH\t: %.1f\tm2/s2\n", (*ptr).efi.esrh);
    printf("\tEBWD\t: %.1f\tm/s\n", (*ptr).efi.ebwd);
    printf("\tSHIP\t:%.1f\n", (*ptr).ship);
    printf("\tSCP\t:%.1f\n", (*ptr).scp);
    printf("Storm motion\n");
    printf("Storm motion vector:(Bunker)\n");
    printf("Left-Mover: ");
    spedir((*ptr).bunker.lstu,(*ptr).bunker.lstv);
    printf("Right-Mover:");
    spedir((*ptr).bunker.rstu,(*ptr).bunker.rstv);
    printf("\n\t\t\tRight-Mover\t\t\t\t\tLeft-Mover\n");
    printf("\t\tSRH0-1\t\tSRH0-2\t\tSRH0-3\t\tSRH0-1\tSRH0-2\tSRH0-3\n");
    printf("\t\tpositv\tnormal\tpositv\tnormal\tpositv\tnormal\t\n");
    printf("Bunker:\t\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",(*ptr).bunker.srh1p,(*ptr).bunker.srh1,(*ptr).bunker.srh2p,(*ptr).bunker.srh2,(*ptr).bunker.srh3p,(*ptr).bunker.srh3,(*ptr).bunker.leftsrh1,(*ptr).bunker.leftsrh2,(*ptr).bunker.leftsrh3);
//    printf("Sharppy:\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",(*ptr).sharppy.srh1p,(*ptr).sharppy.srh1,(*ptr).sharppy.srh2p,(*ptr).sharppy.srh2,(*ptr).sharppy.srh3p,(*ptr).sharppy.srh3,(*ptr).sharppy.leftsrh1,(*ptr).sharppy.leftsrh2,(*ptr).sharppy.leftsrh3);
    printf("\n\tML Wmax Shear (0-6) : %.1f", (*ptr).wxshr);
    printf("\n\tBRN Shear : %.1f", (*ptr).brn_shear);
    printf("\n\tBRN : %.1f", (*ptr).brn_ml);
    printf("\n\tBRN : %.1f", (*ptr).brn_mu);
    printf("\n\tCriticle Angle : %.1f" , (*ptr).bunker.angle);

    return;
}

//Ausgabe in eine HTML-Datei

void htmlausgabe(struct indices *ptr, char *station, struct sounding_data all[], struct tm *t, int j){
    FILE *html;
    int i=0;
    char *name;
    char bild1[16], bild3[15], bild4[14], bild5[15], bild2[16], bild6[16];
    name = (char *) malloc (sizeof(char) * 16);
    /*some string manipulation*/
    memmove(name, "radio",5);
    memmove(name + 5,station,5);   
    memmove(name + 10, ".html\0",6);

    if ( (html = fopen(name,"w+")) == NULL ){
        printf("HTML-Datei konnte nicht geöffnet werden.\n");
        perror(NULL);  
    }

    fprintf(html,"<html>\n  <head>\n    <title>Temp Auswertung</title>\n    <link href=\"https://userpage.fu-berlin.de/mammatus95/temps/dropdown.css\" type=\"text/css\" rel=\"stylesheet\">\n    <link href=\"https://userpage.fu-berlin.de/mammatus95/temps/general.css\" type=\"text/css\" rel=\"stylesheet\">\n    <link href=\"https://userpage.fu-berlin.de/mammatus95/temps/slide.css\" type=\"text/css\" rel=\"stylesheet\">\n");
//CSS-Code:
    fprintf(html,"  </head>\n    <body>\n");
    fprintf(html,"    <h1 class=\"c\">Soundings</h1>    <div class=\"navbar\">      <div class=\"dropdown\">        <button class=\"dropbtn\">Other Analyzes</button>        <div class=\"dropdown-content\">          <a href=\"http://ertel2.uibk.ac.at:8080/raso/\" target=\"_blank\">Ertel</a>          <a href=\"http://weather.uwyo.edu/upperair/\" target=\"_blank\">Upperairmaps&Temps</a>          <a href=\"http://weather.uwyo.edu/upperair/bufrraob.shtml\" target=\"_blank\">Upperairmaps BufrTemps</a>        </div>      </div>      <div class=\"dropdown\">        <button class=\"dropbtn\">Soundings</button>        <div class=\"dropdown-content\">		      <a href=\"radio10393.html\" target=\"Rechtes_Fenster\">Lindenberg</a>          <a href=\"radio11035.html\" target=\"Rechtes_Fenster\">Wien</a>		      <a href=\"radio10184.html\" target=\"Rechtes_Fenster\">Greifswald</a>		      <a href=\"radio10238.html\" target=\"Rechtes_Fenster\">Bergen</a>		      <a href=\"radio10548.html\" target=\"Rechtes_Fenster\">Meiningen</a>		      <a href=\"radio11520.html\" target=\"Rechtes_Fenster\">Praha-Libus</a>		      <a href=\"radio11747.html\" target=\"Rechtes_Fenster\">Prostejov</a>		      <a href=\"radio10410.html\" target=\"Rechtes_Fenster\">Essen</a>        </div>      </div>      <div class=\"dropdown\">        <button class=\"dropbtn\">Forecast Soundings</button>        <div class=\"dropdown-content\">          <a href=\"https://userpage.fu-berlin.de/mammatus95/icon/sounding00.html\">ICON NEST EU</a>                  <a href=\"https://userpage.fu-berlin.de/mammatus95/icon/sounding.html\">ICON NEST EU</a>\n          <a href=\"https://www.tropicaltidbits.com/analysis/models/?model=gfs&region=eu&pkg=mslp_pcpn_frzn&runtime=2018042612&fh=0\">GFS</a>\n          <a href=\"https://apps.ecmwf.int/webapps/opencharts/products/opencharts_vertical-profile-meteogram?&&base_time=202111110000&lat=52.3989&lon=13.0657&station_name=Potsdam&valid_time=202111121200\">IFS</a>    \n        </div>\n      </div>\n");
    fprintf(html,"      <div class=\"dropdown\">\n        <button class=\"dropbtn\">Forecast Soundings</button>        <div class=\"dropdown-content\">\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/icon/sounding00.html\">ICON NEST EU</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/cosmo/soundings.html\">COSMO D2</a>\n        </div>\n      </div>\n");
    fprintf(html,"      <div class=\"dropdown\">\n        <button class=\"dropbtn\">Examples</button>\n        <div class=\"dropdown-content\">          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/example/06.11.2016/radio16245.html\" target=\"Rechtes_Fenster\">Rom 06.11.2016</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/example/05.11.2017/radio16245.html\" target=\"Rechtes_Fenster\">Rom 05.11.2017</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/example/29.10.2018/radio16245.html\" target=\"Rechtes_Fenster\">Rom 29.10.2018</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/example/05.05.2015/radio10184.html\" target=\"Rechtes_Fenster\">Greifswald 05.05.2015 (B&uuml;tzo Tornado)</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/example/09.06.2014/radio10238.html\" target=\"Rechtes_Fenster\">Bergen 09.06.2014 (Pfingsunwetter)</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/example/26.08.2011/radio10238.html\" target=\"Rechtes_Fenster\">Bergen 26.08.2011</a>\n        </div>\n      </div>\n");
    fprintf(html,"      <div class=\"dropdown\">\n        <button class=\"dropbtn\">Klimatologie Lindenberg</button>\n        <div class=\"dropdown-content\">\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/climate/now.html\" target=\"Rechtes_Fenster\">Aktuelles Jahr</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/climate/wind.html\" target=\"Rechtes_Fenster\">Windklimatologie</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/climate/temp.html\" target=\"Rechtes_Fenster\">Temperaturverlauf über das Jahr</a>\n          <a href=\"https://userpage.fu-berlin.de/mammatus95/temps/climate/month.html\" target=\"Rechtes_Fenster\">zeitliche Änderung der Temperatur</a>\n");
    fprintf(html,"\n        </div>\n      </div>\n    </div>\n\n\n<div class=\"m\">\n");
    //fprintf(html,"\n\tMonth: %d Day: %d Time: %d UTC\n Launch Time: %d UTC\n", ((*t).tm_mon + 1), (*t).tm_mday,(*t).tm_hour,(*t).tm_hour - 1);

    //fprintf(html,"\n<h3>Grafiken</h3>\n");
/*Bilder*/
    fprintf(html,"<!-- zuerst das Sichtfenster -->\n<div class=\"cssSlider\">\n\n<!-- die inputs um den Slider zu Steuern -->\n<input name=\"slider\" id=\"slide01\" checked=\"checked\" type=\"radio\">\n<input name=\"slider\" id=\"slide02\" type=\"radio\">\n<input name=\"slider\" id=\"slide03\" type=\"radio\">\n<input name=\"slider\" id=\"slide04\" type=\"radio\">\n<input name=\"slider\" id=\"slide05\" type=\"radio\">\n<input name=\"slider\" id=\"slide06\" type=\"radio\">\n\n<!-- die einzelnen Slides, hier als Liste angelegt -->\n <ul class=\"sliderElements\">\n<li>\n<figure>\n");


    memmove(&bild1[0],station,5);
    memmove(&bild1[0] + 5,"skewT.png\0",10);
    fprintf(html,"\t\t\t<img src=\"%s\"  alt=\"\" width=\"1200\" height=\"400\">\n",&bild1[0]);
    fprintf(html,"\t\t\t<figcaption>Skew T - P log.</figcaption>");
    fprintf(html,"\t\t</figure>\n\t</li>\n\t<li>\n\t\t<figure>");
    memmove(&bild4[0],station,5);   
    memmove(&bild4[0] + 5,"wind.png\0",9);
    fprintf(html,"\t\t\t<img src=\"%s\"  alt=\"\" width=\"1200\" height=\"400\">\n",&bild4[0]);
    fprintf(html,"\t\t\t<figcaption>Hodograph</figcaption>");
    fprintf(html,"\t\t</figure>\n\t</li>\n\t<li>\n\t\t<figure>");
    memmove(&bild2[0],station,5);   
    memmove(&bild2[0] + 5,"nixon.png\0",10);
    fprintf(html,"\t\t\t<img src=\"%s\"  alt=\"\" width=\"1200\" height=\"400\">\n",&bild2[0]);
    fprintf(html,"\t\t\t<figcaption>storm relative Hodograph</figcaption>");
    fprintf(html,"\t\t</figure>\n\t</li>\n\t<li>\n\t\t<figure>");
    memmove(&bild5[0],station,5);
    memmove(&bild5[0] + 5, "thermo.png\0",11);
    fprintf(html,"\t\t\t<img src=\"%s\"  alt=\"\" width=\"1200\" height=\"400\">\n",&bild5[0]);
    fprintf(html,"\t\t\t<figcaption>Stuve Diagramm: Temperatur (Blau)- und Taupunktskurve (Gr&uuml;n), sowie das MU-Parcel(Rot) und EFI (blaue Horizontale Linien)</figcaption>");
    fprintf(html,"\t\t</figure>\n\t</li>\n\t<li>\n\t\t<figure>");
    memmove(&bild3[0],station,5);   
    memmove(&bild3[0] + 5,"theta.png\0",10); 
    fprintf(html,"\t\t\t<img src=\"%s\"  alt=\"\" width=\"1200\" height=\"400\">\n",&bild3[0]);
    fprintf(html,"\t\t\t<figcaption>Verlauf Potentieller und Pseudopotentieller Temperatur mit dem Druck als Vertikalkoordinate.</figcaption>");
    fprintf(html,"\t\t</figure>\n\t</li>\n\t<li>\n\t\t<figure>");
    memmove(&bild6[0],station,5);   
    memmove(&bild6[0] + 5,"ri.png\0",7);
    fprintf(html,"\t\t\t<img src=\"%s\"  alt=\"\" width=\"1200\" height=\"400\">\n",&bild6[0]);
    fprintf(html,"\t\t\t<figcaption></figcaption>");
    
    fprintf(html,"\t\t</figure>\n\t</li>\n</ul>\n\t<!-- Eine Steuerung -->\n\t<ul class=\"sliderControls\">\n\t\t<li><label for=\"slide01\">1</label></li>\n\t\t<li><label for=\"slide02\">2</label></li>\n\t\t<li><label for=\"slide03\">3</label></li>\n\t\t<li><label for=\"slide04\">4</label></li>\n\t\t<li><label for=\"slide05\">5</label></li>\n\t\t<li><label for=\"slide06\">6</label></li>\n\t</ul>\n</div>\n\n");

    //memmove(&bild5[0],station,5);   
    //memmove(&bild5[0] + 5,"ver.png\0",8);

    fprintf(html,"\t</br></br><h3>Thermodynamsiche Indizes</h3></br>\n");
    fprintf(html,"\tLapse Rate zwischen 850 & 500\t\t:%.1f\tK/km</br>\n", (*ptr).lapse1);
    fprintf(html,"\tLapse Rate zwischen 700 & 500\t\t:%.1f\tK/km</br>\n", (*ptr).lapse2);
    fprintf(html,"\tLapse Rate zwischen 925 & 850\t\t:%.1f\tK/km</br>\n", (*ptr).lapse3);  
    fprintf(html,"\tKO-Index:\t:%.1f\tK\t <a href=\"http://www1.wetter3.de/ko_index.html\" target=\"_blank\" ><font size=\"2px\" color=\"#FF0000\">Erkl&auml;rung KO-Index</font></a> </br>\n", (*ptr).koindex);
    fprintf(html,"\tK-Index:\t:%.1f\tK\t<a href=\"http://glossary.ametsoc.org/wiki/Stability_index\" target=\"_blank\" ><font size=\"2px\" color=\"#FF0000\">Erkl&auml;rung Stabilit&auml;ts Indizes</font></a> \n\t</br>\n", (*ptr).kindex);
    fprintf(html,"\tTT-Index:\t:%.1f\tK\t Miller (1972)\n \t</br>\n", (*ptr).tt);
    fprintf(html,"\tTemperatur Maximum\t:%.1f\t&ordm;C\n", (*ptr).tx);
    fprintf(html,"\tMean Mix Ratio\t:%.1f\tg/kg\n\n", (*ptr).mmix);
    fprintf(html,"\t<h4>CAPE</h4>\n");
    fprintf(html,"\t<table>\n\t\t<tr>\n\t\t\t<th></th><th>CAPE</th><th>CIN</th><th>LCL High</th><th>LFC High</th><th>EQ High</th>\n\t\t</tr>\n");
    fprintf(html,"\t\t<tr>\n\t\t\t<th></th><th>J/kg</th><th>J/kg</th><th>m</th><th>m</th><th>m</th>\n\t\t</tr>\n");
    fprintf(html,"\t\t<tr>\n\t\t\t<td align=\"center\">Surface base Parcel</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td>\n\t\t</tr>\n",(*ptr).sb.cape,(*ptr).sb.cin,(*ptr).sb.lclh,(*ptr).sb.lfch,(*ptr).sb.eqh);
    fprintf(html,"\t\t<tr>\n\t\t\t<td align=\"center\">Most Unstable Parcel</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td>\n\t\t</tr>\n",(*ptr).mu.cape,(*ptr).mu.cin,(*ptr).mu.lclh,(*ptr).mu.lfch,(*ptr).mu.eqh);
    fprintf(html,"\t\t<tr>\n\t\t\t<td align=\"center\">Mixed Layer Parcel</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td><td align=\"center\">%.1f</td>\n\t\t</tr>\n",(*ptr).ml.cape,(*ptr).ml.cin,(*ptr).ml.lclh,(*ptr).ml.lfch,(*ptr).ml.eqh);
    fprintf(html,"\n\t\t</table>\n</br></br>");

    fprintf(html,"\n\n\t<h5>Wind</h5>\n");
    fprintf(html,"\tbulk Richardson number\t:%.1f\t <a href=\"http://tornado.sfsu.edu/geosciences/classes/m201/buoyancy/SkewTMastery/mesoprim/skewt/brn.htm\" target=\"_blank\" ><font size=\"2px\" color=\"#FF0000\">    Erkl&auml;rung BRN</font></a></br>\n",(*ptr).brn_ml);
    fprintf(html,"\tBRN Shear \t:%.1f\tm2/s2       Values of 35-40 m2/s2 or greater have been associated with supercells.</br>\n",(*ptr).brn_shear);
    fprintf(html,"\tBulk Shear:</br>\n");
    fprintf(html,"\tShear 0-500 m\t:%.1f\tm/s</br>\n",(*ptr).shear500m);
    fprintf(html,"\tShear 0-1 km\t:%.1f\tm/s</br>\n",(*ptr).shear1);
    fprintf(html,"\tShear 0-3 km\t:%.1f\tm/s</br>\n",(*ptr).shear3);
    fprintf(html,"\tShear 0-6 km\t:%.1f\tm/s</br>\n",(*ptr).shear6);
    fprintf(html,"\tShear 0-8 km\t:%.1f\tm/s\t(>25 m/s long lived)</br>\n",(*ptr).shear8);
    fprintf(html,"\tHelicity 0-500 m\t:%.1f\tm2/s2</br>\n", (*ptr).hel500m);
    fprintf(html,"\tHelicity 0-1 km\t:%.1f\tm2/s2</br>\n", (*ptr).hel1);
    fprintf(html,"\tHelicity 0-2 km\t:%.1f\tm2/s2</br>\n", (*ptr).hel2);
    fprintf(html,"\tHelicity 0-3 km\t:%.1f\tm2/s2</br>\n", (*ptr).hel3);
    fprintf(html,"\tSRH 0-1 km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.srh1);
    fprintf(html,"\tSRH 0-2 km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.srh2);
    fprintf(html,"\tSRH 0-3 km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.srh3);

    fprintf(html,"</br>Effektive Inflow Layer</br>\n");
    fprintf(html,"\tESRH\t: %.1f\tm2/s2</br>\n", (*ptr).efi.esrh);
    fprintf(html,"\tEBWD\t: %.1f\tm/s \t <a href=\"http://www.spc.noaa.gov/exper/mesoanalysis/help/begin.html\" target=\"_blank\" ><font size=\"2px\" color=\"#FF0000\">Effective Bulk Wind Difference</font></a> </br>\n", (*ptr).efi.ebwd);
    fprintf(html,"\tSCP\t:%.1f\t</br>\n", (*ptr).scp);
    fprintf(html,"\tSHIP\t:%.1f\t</br>\n", (*ptr).ship);
    fprintf(html,"\tWet Bulb Zero\t:%.1f m\t %.1f ft\t (7000ft and 10.500ft) <a href=\"https://forecast.weather.gov/glossary.php?word=wet%%20bulb%%20zero\" target=\"_blank\" ><font size=\"2px\" color=\"#FF0000\">WBZ</font></a>\t</br>\n", (*ptr).wetbulbzero,3.28084*(*ptr).wetbulbzero);
    fprintf(html,"\tML Wmax Shear\t:%.1f\tm2/s2</br>\n", (*ptr).wxshr);
    fprintf(html,"\tEnergy-Helicity Index\t:%.1f\tm2/s2\t<a href=\"http://www.spc.noaa.gov/exper/mesoanalysis/help/help_ehi3.html\" target=\"_blank\" ><font size=\"2px\" color=\"#FF0000\"> EHI values greater than 1-2 have been associated with significant tornadoes in supercells.</font></a></br>\n", (*ptr).ehi);

    fprintf(html,"\t</br></br><h3>Hauptdruckfl&auml;chen</h3></br>\n");
    fprintf(html,"\t<table>\n\t\t<tr>\t\t\t<th>Druck</th>\n\t\t\t<th>Temperatur</th>\n\t\t\t<th>Taupunkt</th>\n\t\t\t<th>Windgeschwindigkeit</th>\n\t\t\t<th>Windrichtung</th>\n\t\t</tr>\n");
    while (i < j){
        if (all[i].druck == 850.0)
            fprintf(html,"\t\t<tr>\n\t\t\t<td align=\"center\">%d</td>\n\t\t\t<td align=\"center\">%.1f C</td>\n\t\t\t<td align=\"center\">%.1f C</td>\n\t\t\t<td align=\"center\">%.1f kn</td>\n\t\t\t<td align=\"center\">%.1f</td>\n\t\t</tr>\n",(int) all[i].druck,all[i].temp,all[i].dewpoint,all[i].wspd,all[i].wdir);
        if (all[i].druck == 700.0)
            fprintf(html,"\t\t<tr>\n\t\t\t<td align=\"center\">%d</td>\n\t\t\t<td align=\"center\">%.1f C</td>\n\t\t\t<td align=\"center\">%.1f C</td>\n\t\t\t<td align=\"center\">%.1f kn</td>\n\t\t\t<td align=\"center\">%.1f</td>\n\t\t</tr>\n",(int) all[i].druck,all[i].temp,all[i].dewpoint,all[i].wspd,all[i].wdir);
        if (all[i].druck == 500.0)
            fprintf(html,"\t\t<tr>\n\t\t\t<td align=\"center\">%d</td>\n\t\t\t<td align=\"center\">%.1f C</td>\n\t\t\t<td align=\"center\">%.1f C</td>\n\t\t\t<td align=\"center\">%.1f kn</td>\n\t\t\t<td align=\"center\">%.1f</td>\n\t\t</tr>\n",(int) all[i].druck,all[i].temp,all[i].dewpoint,all[i].wspd,all[i].wdir);
        if (all[i].druck == 300.0)
            fprintf(html,"\t\t<tr>\n\t\t\t<td align=\"center\">%d</td>\n\t\t\t<td align=\"center\">%.1f C</td>\n\t\t\t<td align=\"center\">%.1f C</td>\n\t\t\t<td align=\"center\">%.1f kn</td>\n\t\t\t<td align=\"center\">%.1f</td>\n\t\t</tr>\n",(int) all[i].druck,all[i].temp,all[i].dewpoint,all[i].wspd,all[i].wdir);
        i++;
    }
    fprintf(html,"\t<table>\n");

    fprintf(html,"\t<h2>Erg&auml;nzungen</h2>\n");
    fprintf(html,"\tbulk Richardson number MU Parcel\t:%.1f\t</br>\n",(*ptr).brn_mu);
    fprintf(html,"\tSRH positiv 0-1km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.srh1p);
    fprintf(html,"\tSRH positiv 0-2km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.srh2p);
    fprintf(html,"\tSRH positiv 0-3km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.srh3p);
    fprintf(html,"\tSRH left 0-1km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.leftsrh1);
    fprintf(html,"\tSRH left 0-2km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.leftsrh2);
    fprintf(html,"\tSRH left 0-3km\t:%.1f\tm2/s2</br>\n", (*ptr).bunker.leftsrh3);
    fprintf(html,"\tStreamwise 0-1km\t:%.1f\t1/s</br>\n", (*ptr).bunker.stw1);
    fprintf(html,"\tStreamwise 0-2km\t:%.1f\t1/s</br>\n", (*ptr).bunker.stw2);
    fprintf(html,"\tStreamwise 0-3km\t:%.1f\t1/s</br>\n", (*ptr).bunker.stw3);
    fprintf(html,"\tStrom Relative Wind Shear 0-1 km\t:%.1f\tm/s</br>\n", (*ptr).bunker.sr_shear1);
    fprintf(html,"\tStrom Relative Wind Shear 0-3 km\t:%.1f\tm/s</br>\n", (*ptr).bunker.sr_shear3);
    fprintf(html,"\tStrom Relative Wind Shear 0-6 km\t:%.1f\tm/s</br>\n", (*ptr).bunker.sr_shear6);
    fprintf(html,"\tStrom Relative Wind Shear 2-4 km\t:%.1f\tm/s</br>\n", (*ptr).bunker.sr_shear24);
    fprintf(html,"\tLength of the Hodograph between 0-1 km\t:%.1f\tm/s</br>\n", (*ptr).length1);
    fprintf(html,"\tLength of the Hodograph between 0-3 km\t:%.1f\tm/s</br>\n", (*ptr).length3);
    fprintf(html,"\tLength of the Hodograph between 0-6 km\t:%.1f\tm/s</br>\n", (*ptr).length6);
    fprintf(html,"\tLength of the Hodograph between 0-8 km\t:%.1f\tm/s</br>\n", (*ptr).length8);
    fprintf(html,"\tLength of the Hodograph between 2-8 km\t:%.1f\tm/s</br>\n", (*ptr).length28);


    fprintf(html,"    </div>\n  </body>\n</html>");
    fclose(html);
}

#endif /*AUSGABE_H*/
