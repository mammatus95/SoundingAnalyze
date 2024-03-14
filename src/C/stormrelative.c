/* Wyoming Sounding as input with strom motion as argv*/


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
        //for(i=0;i < n; i++)
        //    printf("%d : %p : %s\n",i,token[i],token[i]);
        if ( n > 10){
            (*struct_ptr).druck = atoi(token [0]);//printf("%d",(*struct_ptr).druck);
            (*struct_ptr).high = atoi(token [1]);
            (*struct_ptr).temp = strtod(token [2],NULL);
            (*struct_ptr).dewpoint = strtod(token [3],NULL);
            (*struct_ptr).wdir = strtod(token [6],NULL);
            (*struct_ptr).wspd = KT2MS(strtod(token [7],NULL));
        }
    } else {return -1;}
    free (buffer);
    return 0;
}

int main (int argc, char *argv[]) {

    FILE *quelle, *data;
    int j = 0, i = 0, anzahl=200, ret=0;
    char buffer[90], aufruf[120] = {"python3 plottemp.py "}, station[6];

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
            printf("Anzahl : %d",anzahl);
        } else {
            exit(0);
        }
    } else {
        exit(0);
    }

    struct sounding_data drk[anzahl];


    if( (data=fopen("radiodata.txt","w+")) == NULL) {
        fprintf(stderr, "Kann radiodata.txt nicht oeffnen\n");
        return EXIT_FAILURE;
    }
    /*Beginn der Auswertung der Radiosondendaten*/

    for(j=0;j<=anzahl;j++){
        if ( 0 == (zerlegen(quelle,&drk[j])) ){
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
        }else{
            break;
        }
    }
    anzahl=j;

    severe.ml = *(mixed_layer(drk,&severe));
    ret = parcel_wrapper (drk,&severe.mu,&severe.sb);

    /*interesting*/
    for(int k=0;k<=anzahl;k++){
        if (drk[k].cape == 0)
            drk[k].cin = 0;
    }

