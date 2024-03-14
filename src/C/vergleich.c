#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct info{
    int druck;
    int high;
    float temp;
    float dewpoint;
    float pottemp;
    float wdir;
    float wspd;
};

static char *token[30];

/*Einlesen der Input-File*/

/*Von Wyoming*/
int zerlegen(FILE *quelle, struct info *struct_ptr){
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
        if ( n > 3){
            (*struct_ptr).druck = atoi(token [0]);//printf("%d",(*struct_ptr).druck);
            (*struct_ptr).high = atoi(token [1]);
            (*struct_ptr).temp = strtod(token [2],NULL);
            (*struct_ptr).dewpoint = strtod(token [3],NULL);
            (*struct_ptr).wdir = strtod(token [6],NULL);
            (*struct_ptr).wspd = strtod(token [7],NULL);
        }
    } else {return -1;}
    free (buffer);
    return 0;
}

//Ausgabe in die Datei vergleich.txt
void fileoutput (struct info *struct_ptr, struct info *ptr2, FILE *data){

    fprintf(data,"%d,", (*struct_ptr).high);
    fprintf(data,"%d,", (*struct_ptr).druck);
    fprintf(data,"%.1f,", (*struct_ptr).temp);
    fprintf(data,"%.1f,",(*struct_ptr).dewpoint);
    fprintf(data,"%d,", (*ptr2).high);
    fprintf(data,"%d,", (*ptr2).druck);
    fprintf(data,"%.1f,", (*ptr2).temp);
    fprintf(data,"%.1f\n",(*ptr2).dewpoint);

    return;
}

int main (int argc, char *argv[]) {

    FILE *ziel, *data1, *data2;
    int j = 0, anzahl=200, ret=0;


    struct info temp1[anzahl], temp2[anzahl];


    /*Dteien öffnen*/

    if ( argc == 3 ){
        if ( ( strlen(argv[1]) < 13 ) && ( strlen(argv[2]) < 13 ) ){
            if ( (data1 = fopen(argv[1],"r")) == NULL ){
                printf("Datei konnte nicht geöffnet werden.\n");
                perror(NULL);
                exit(0); 
            }
            if ( (data2 = fopen(argv[2],"r")) == NULL ){
                printf("Datei konnte nicht geöffnet werden.\n");
                perror(NULL);
                exit(0); 
            }
        } else {
            exit(0);
        }
        
    } else {
        exit(0);
    }

    if( (ziel=fopen("vergleich.txt","w+")) == NULL) {
        fprintf(stderr, "Kann radiodata.txt nicht oeffnen\n");
        return EXIT_FAILURE;
    }

    while(j<=anzahl){
        if( (temp1[j].druck > 100 ) || (temp2[j].druck > 100) ){
            zerlegen(data1,&temp1[j]);
            zerlegen(data2,&temp2[j]);
            fileoutput(&temp1[j],&temp2[j],ziel);
        } else {
            break;
        }
    }

    fclose(data1);
    fclose(data2);
    fclose(ziel);

    //system("python2 vergleich.py");

    return 0;
}
