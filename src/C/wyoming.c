#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sounding.h"
#include <unistd.h>
#include <dirent.h>



static char *token[30];

/*Einlesen der Input-File*/

/*Von Wyoming*/
int zerlegen(FILE *quelle, struct sounding_data *struct_ptr){
    char *buffer = (char *) malloc (300);
    int n =0;

    if( (fgets(buffer, 298, quelle)) != NULL){
        n = 0;
        token[n] = strtok(buffer," ");
        while(token[n] != NULL){
            n++;
            token[n] = strtok(NULL," ");
        }
        //if ( n < 4){
        //for(int i=0;i < n; i++)
        //    printf("%d : %p : %s\n",i,token[i],token[i]);
        //}
        if (n > 3){
            (*struct_ptr).druck = strtod(token [0],NULL);//printf("%d",(*struct_ptr).druck);
            (*struct_ptr).high = strtod(token [1],NULL);
            (*struct_ptr).temp = strtod(token [2],NULL);
            (*struct_ptr).dewpoint = strtod(token [3],NULL);
            (*struct_ptr).wdir = strtod(token [6],NULL);
            (*struct_ptr).wspd = strtod(token [7],NULL); //KT2MS(strtod(token [7],NULL));
            free (buffer);return 0;
        } else {
            free (buffer);return 1;
        }
    } else {free (buffer);return -1;}
    free (buffer);
    return 0;
}
/*Turm*/
//500 Elemte ca.
/*
const char  such[] = { "severe 10393" };// 07.06.17 12 UTC Lindenberg" };

int zerlegen(FILE *quelle, struct sounding_data *struct_ptr){
    char *buffer = (char *) malloc (300);
    int n,i;
    if( (fgets(buffer, 298, quelle)) != NULL){
        if ( (strstr(buffer, "2222")) != NULL)
            return 2;

        n = 0;
        token[n] = strtok(buffer," ");
        while(token[n] != NULL){
            n++;
            token[n] = strtok(NULL," ");
        }
        for(i=0;i < n; i++)
            printf("%d : %p : %s\n",i,token[i],token[i]);
        if ( n > 3){
            (*struct_ptr).druck = atoi(token [0]);//printf("%d",(*struct_ptr).druck);
            (*struct_ptr).high = atoi(token [1]);
            (*struct_ptr).severe = strtod(token [2],NULL);
            (*struct_ptr).dewpoint = strtod(token [3],NULL);
        }
    } else {return -1;}
    free (buffer);
    return 0;
}

int zerlegenwind (FILE *quelle, struct sounding_data *struct_ptr){
    char *buffer = (char *) malloc (300);
    int n,i;
    if( (fgets(buffer, 298, quelle)) != NULL){
        if ( (strstr(buffer, "9999")) != NULL)
            return 2;

        n = 0;
        token[n] = strtok(buffer," ");
        while(token[n] != NULL){
            n++;
            token[n] = strtok(NULL," ");
        }
        for(i=0;i < n; i++)
            printf("%d : %p : %s\n",i,token[i],token[i]);
        if ( n > 2){
            //(*struct_ptr).druck = atoi(token [0]);//printf("%d",(*struct_ptr).druck);
            (*struct_ptr).wdir = strtod(token [1],NULL);
            (*struct_ptr).wspd = strtod(token [2],NULL);

        }
    } else {return -1;}
    free (buffer);
    return 0;
}
*/

int main (int argc, char *argv[]) {

    FILE *quelle, *data;
    int j = 0, i = 0, anzahl=200, ret=0, ret_mu=0;
    char buffer[90], aufruf[164] = {"/home/mammatus95/Dokumente/miniconda3/bin/python3 plottemp.py "}, station[6];

    struct indices severe;
    struct tm time;
    
    //printf("%p\n\n", &drk[0]);

    /*Dteien öffnen*/

    if ( argc == 2 ){
        if ( strlen(argv[1]) < 10 ){
            if ( (quelle = fopen(argv[1],"r")) == NULL ){
                printf("Datei konnte nicht geöffnet werden.\n");
                perror(NULL);
                exit(0); 
            }
            strncpy(&station[0],argv[1],5);
            strncpy(&station[0]+5,"\0",1);
        } else {
            exit(0);
        }
    } else if (argc == 6){
        anzahl = atoi(argv[1]);
        if ( strlen(argv[2]) < 10 ){
            if ( (quelle = fopen(argv[2],"r")) == NULL ){
                printf("Datei konnte nicht geöffnet werden.\n");
                perror(NULL);
                exit(0);
            }
            strncpy(&station[0],argv[2],5);
            strncpy(&station[0]+5,"\0",1);
            /*time*/
            time.tm_mon = ( atoi(argv[3]) - 1);
            time.tm_mday = atoi(argv[4]);
            time.tm_hour = atoi(argv[5]);
            //printf("Anzahl : %d\n",anzahl);
        } else {
            exit(0);
        }
    } else {
        exit(0);
    }

    struct sounding_data drk[anzahl];


/*
    struct sounding_data drk[500];
    if ( (quelle = fopen("20170607_12.TLP", "r")) == NULL ){
        printf("Datei konnte nicht geöffnet werden.\n");
        perror(NULL);  
    }
*/   
    if( (data=fopen("radiodata.txt","w+")) == NULL) {
        fprintf(stderr, "Kann radiodata.txt nicht oeffnen\n");
        return EXIT_FAILURE;
    }
    /*Beginn der Auswertung der Radiosondendaten*/
    j=0;
    while(j<=anzahl){
        ret = zerlegen(quelle,&drk[j]);
        //printf("%d\t%d\n",j,ret);
        if ( ret == 0 ){
            calculate(&drk[j]);
            if (j < 1){
                drk[0].lapse = NAN;
                drk[0].ushr = 0;
                drk[0].vshr = 0;
                drk[0].cape = 0;
                drk[0].cin = 0;
                drk[0].agl = 0;
            } else {
                drk[j].agl=drk[j].high-drk[0].high;
                windshear(&drk[j-1],&drk[j]);
                lapserate(&drk[j-1],&drk[j]);
            }
            drk[j].cape = 0;
            drk[j].cin = 0;
            if (drk[j].druck < 50.0 ) {
                break;
            }
            j++;
        } else if (ret == 1 ) {
            continue;
        }else{
            break;
        }

    }
    anzahl=j;ret=0;
    //printf("\n%d\n",anzahl);

    /*Calculation of Must Unstable CAPE*/
    //severe.mu = (*cape(drk));
    //ret = severe.mu.index;
    //severe.ml = (*ml_cape(drk));

    severe.ml = *(mixed_layer(drk,&severe));
    ret = parcel_wrapper (drk,&severe.mu,&severe.sb); ret_mu = ret;
    /*interesting*///? 
    for(int k=0;k<=anzahl;k++){
        if( (ret==0) && (drk[k].wetbulb <= 0.0) ){
            severe.wetbulbzero = drk[k].agl;
            ret=1;
        }
        if (drk[k].cape == 0)
            drk[k].cin = 0;
    }
    ret=0;
    /*Calculation of Helicity*/

    severe.hel500m = helicity(drk,500);
    severe.hel1 = helicity(drk,1000);
    severe.hel2 = helicity(drk,2000);
    severe.hel3 = helicity(drk,3000);

    /*Calculation the shear*/

    severe.shear500m = scherung(drk,500);
    severe.shear1 = scherung(drk,1000);
    severe.shear3 = scherung(drk,3000);
    severe.shear6 = scherung(drk,6000);
    severe.shear8 = scherung(drk,8000);

    severe.length1 = length(drk,0,1000);
    severe.length3 = length(drk,0,3000);
    severe.length6 = length(drk,0,6000);
    severe.length8 = length(drk,0,8000);
    severe.length28 = length(drk,2000,8000);

    //Berechnung der Windauswertung

    

    severe.bunker = (*parcel_bunkers_motion(drk,&severe.mu,j));
    for (i=0; i <= j;i++){
        stormrelativwind(&severe.bunker, &drk[i]);
    }

    severe.efi = (*effective_inflow_layer(drk,severe.mu.eqh,j));

    severe.bunker.srh1 = srh(drk,1000,1);
    severe.bunker.srh2 = srh(drk,2000,1);
    severe.bunker.srh3 = srh(drk,3000,1);

    severe.bunker.srh1p = srh(drk,1000,3);
    severe.bunker.srh2p = srh(drk,2000,3);
    severe.bunker.srh3p = srh(drk,3000,3);

    severe.bunker.leftsrh1 = srh(drk,1000,2);
    severe.bunker.leftsrh2 = srh(drk,2000,2);
    severe.bunker.leftsrh3 = srh(drk,3000,2);    

    severe.bunker.stw1 = streamwise(drk,1000);
    severe.bunker.stw2 = streamwise(drk,2000);
    severe.bunker.stw3 = streamwise(drk,3000);

    severe.bunker.sr_shear1 = stormrelativ_shear(drk,1000);
    severe.bunker.sr_shear3 = stormrelativ_shear(drk,3000);
    severe.bunker.sr_shear6 = stormrelativ_shear(drk,6000);
    severe.bunker.sr_shear24 = stormrelativ_shear2(drk,4000,2000);

    /*Errechung weiterer Indizes*/
    research(drk,&severe,anzahl);

    severe.ehi = ehi(severe.bunker.srh3p, severe.mu.cape);

    severe.ship = ship (drk, severe.shear6, ret, &severe.mu.cape, severe.lapse2,j);

    severe.scp = scp(&severe.mu.cape, severe.efi.esrh, severe.efi.ebwd);

    brnshear(drk,&severe);

    severe.wxshr = wxshear(&severe.ml.cape,&severe.shear6);

    severe.tx = t_max(drk);

    //Ausgabe in die Datei
    for(i=0;i<anzahl;i++){
        fileoutput(&drk[i],data);
        //ausgabe1(&drk[i],i);
    }
    fclose(data);

    //Ausgabe Indizes
    //ausgabe2(&severe,&time);
    htmlausgabe(&severe,&station[0],drk,&time,j);

    /*Vorbereitung Kommandozeilenargument für Python-Skript*/

    if ( snprintf(buffer, sizeof buffer, "%s %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f", &station[0], ret_mu, severe.mu.cape, severe.mu.lclp, severe.mu.lclt, severe.bunker.rstu, severe.bunker.rstv, severe.bunker.lstu, severe.bunker.lstv, severe.efi.pbot, severe.efi.ptop, severe.wxshr, severe.brn_ml, severe.lapse1, severe.mmix)< 0) 
        printf("Hat nicht funktioniert.");
    
    strncat(aufruf,buffer,strlen(buffer));
    
    //puts(buffer);
    //puts(aufruf);
    //printf("%lu\n",strlen(aufruf));
    //printf("%d",ret);

    system(aufruf);

    return 0;
}
