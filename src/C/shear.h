#ifndef SHEAR_H
#define SHEAR_H


/*Function to analyse the wind profil of a sounding*/



void uvwind (struct sounding_data *struct_ptr){
    /* 
    This function convert wind speed and direction in to u,v component.
    The Winddirection is given in degree, 
    */
    double val;
    val = PI / 180.0; //degrees are converted to radians
    (*struct_ptr).u = (-1) * (*struct_ptr).wspd * sin((*struct_ptr).wdir * val);
    (*struct_ptr).v = (-1) * (*struct_ptr).wspd * cos((*struct_ptr).wdir * val);
    
    return;
}

void spedir (float u, float v){
    /*
    Convert the u,v component's back to speed and direction.
    */
    double val, direction, speed;
    val = 180 / PI;
    if ( (direction = atan2(-u,-v) * val) < 0) 
        direction += 360;
    speed = sqrt( SQR(u) + SQR(v) );
    printf("\tRichtung:%.1lf\tGeschwindigkeit:%.1lf\n",direction,speed);
    return;
}

void windshear (struct sounding_data *bot, struct sounding_data *top){
    /*
    This function gives back the wind difference between two layers.
    This equal with the thermal wind.
    Parameters are the u,v Komponente of the top and bottom layer.
    */
    (*top).ushr = (*top).u - (*bot).u;
    (*top).vshr = (*top).v - (*bot).v;
    return;
}

float shear (float ubot, float vbot, float utop, float vtop){
    /*
    This function returns the magnitut of the wind difference between two layers.
    Parameters are the u,v Komponente of the top and bottom layer.
    */
    float du, dv;
    du = utop - ubot;
    dv = vtop - vbot;
    return sqrt( SQR(du) + SQR(dv) );
}

float scherung (struct sounding_data drk[], float height){
    /*Bereite shear() vor.*/
    int i=0;
    float dz=0,u_top=0,v_top=0;
    while (drk[i].agl < height){
        i++;
    }
    //wmean
    dz = drk[i].agl - drk[i+1].agl; 
    u_top = drk[i].u * ((height - drk[i+1].agl)/dz) +  drk[i+1].u * ((drk[i].agl - height)/dz);
    v_top = drk[i].v * ((height - drk[i+1].agl)/dz) +  drk[i+1].v * ((drk[i].agl - height)/dz);

    return shear(drk[0].u, drk[0].v, u_top, v_top);
}

/*Storm relativ wind shear*/
float stormrelativ_shear (struct sounding_data drk[], float height){

    int i=0;
    float dz=0,sru_top=0,srv_top=0;
    while (drk[i].agl < height){
        i++;
    }
    //wmean
    dz = drk[i].agl - drk[i+1].agl;
    sru_top = drk[i].srur * ((height - drk[i+1].agl)/dz) +  drk[i+1].srur * ((drk[i].agl - height)/dz);
    srv_top = drk[i].srvr * ((height - drk[i+1].agl)/dz) +  drk[i+1].srvr * ((drk[i].agl - height)/dz);

    return shear(drk[0].srur, drk[0].srvr, sru_top, srv_top);
}

float stormrelativ_shear2 (struct sounding_data drk[], float top, float bot){

    int i=0;
    float dz=0,sru_top=0,srv_top=0,sru_bot=0,srv_bot=0;

    while (drk[i].agl < bot){
        i++;
    }

    //wmean
    dz = drk[i].agl - drk[i+1].agl;
    sru_bot = drk[i].srur * ((bot - drk[i+1].agl)/dz) +  drk[i+1].srur * ((drk[i].agl - bot)/dz);
    srv_bot = drk[i].srvr * ((bot - drk[i+1].agl)/dz) +  drk[i+1].srvr * ((drk[i].agl - bot)/dz);

    while (drk[i].agl < top){
        i++;
    }
    //wmean
    dz = drk[i].agl - drk[i+1].agl;
    sru_top = drk[i].srur * ((top - drk[i+1].agl)/dz) +  drk[i+1].srur * ((drk[i].agl - top)/dz);
    srv_top = drk[i].srvr * ((top - drk[i+1].agl)/dz) +  drk[i+1].srvr * ((drk[i].agl - top)/dz);

    return shear(sru_bot,srv_bot, sru_top, srv_top);
}



/* Calculates a pressure-weighted mean wind through a layer. */
/*
schoner fehler
float mean_wind_u (struct sounding_data drk[], float bot_high, float top_high){
    int i=0;
    float sum=0 , sum_ps=0;
    while ( (drk[i].agl >= bot_high) && (drk[i].agl <= top_high) ){
        sum += drk[i].u * drk[i].druck;
        sum_ps += drk[i].druck;
        i++;
    }
    printf("\tMean: %.1f\t%.1f\n",sum,sum_ps);
    return sum/sum_ps;
}
*/

float mean_wind_u (struct sounding_data drk[], float bot_high, float top_high){
    int i=0;
    float sum=0 , sum_ps=0;
    while ( drk[i].agl <= top_high ){
        if (drk[i].agl >= bot_high){
            sum += drk[i].u * drk[i].druck;
            sum_ps += drk[i].druck;
        }
        i++;
    }
    return sum/sum_ps;
}


float mean_wind_v (struct sounding_data drk[], float bot_high, float top_high){
    int i=0;
    float sum=0 , sum_ps=0;
    while (drk[i].agl <= top_high){
        if (drk[i].agl >= bot_high){ 
            sum += drk[i].v * drk[i].druck;
            sum_ps += drk[i].druck;
        }
        i++;
    }
    return sum/sum_ps;
}

float helicity ( struct sounding_data all[], float high){
    /*
    The Helicity H is given by (- (dv)/(dz))*u + ((du)/(dz))v for infesitimal difference between two layers
    We approximate this multiply by the layer difference dz:
    H =  sum of all layer ( du * v - dv * u )
    For the SRH become:
    SRH = sum of all layer ( du * srv - dv * sru )
    
    Paramter/Varibales
    dv = vshr = difference between the v-component of two layers / shear vectors of each layer
    du = ushr = difference between the u-component of two layers / shear vectors of each layer
    srv = v-component Storm relativ wind
    sru = u-component Storm relativ wind

    Unit of H is m2 / s2, so we have convert kn in m/s

    */

    float h=0, nH=0, H=0;
    int i=0;

    while (all[i].agl < high){

        //h = ( (all[i].ushr * all[i].v ) - ( all[i].vshr * all[i].u ) );
        h = all[i].ushr * ((all[i].v + all[i+1].v)/ (float) 2) - all[i].vshr * ((all[i].u + all[i+1].u)/ (float) 2);
        //printf("Helicity: %f : %d\n", h, all[i].druck);

        if(h<0){
            nH += h;
        }else if (h >= 0){
            H += h;
        }    
        i++;
        h=0;
    }
    //wmean
    /*Eigentlich musste hier und bei allen weiteren Berechnung noch die ein weight mean fur mehr accurateness angewendet werden.*/
    //printf("negativ: %f\tpositiv: %f\n", nH, H);
    return H + nH;
}

struct stormmotion *parcel_bunkers_motion (struct sounding_data all[], struct parcel *mu, int anzahl){
    /*with struct sounding_data all[] we have whole struct of the "drk" array.*/
    /*The Bunkers motion provides results similar to the "30 degrees right and 75% of the mean wind speed" estimates for typical southwest flow regimes.
    Graphically, the Bunkers storm motion can be estimated by 1) plotting the shear vector from the 0-500 m AGL mean wind to the 5500-6000 m AGL mean wind
    on a hodograph, 2) plotting the 0-6 km mean wind (pressure weighted), and 3) drawing a vector (of 7.5 m/s magnitude) perpendicular to the shear vector
    from the 0-6 km mean wind. A perpendicular vector to the right represents the right (cyclonic in northern hemisphere) supercell motion, and a left
    vector represents the left (anticyclonic in northern hemisphere) supercell motion.
    */
    /*Bunkers Storm Motion Vectors*/

/*
    mucape = (*mu).cape;
    mucinh = (*mu).cin;
    muel = (*mu).eqh;
    base = drk[0].high;
*/
    printf("\n\nMessages from Bunker:\n");
    struct stormmotion *bu;
    bu = malloc(sizeof(struct stormmotion));
    float d = 7.5, tmp=0; //Deviation value emperically derived at 7.5 m/s
    float shru = NAN, shrv = NAN, meanu6 = NAN, meanv6 = NAN;
    //float dz=0,u_top=0,v_top=0;

    float depth = NAN, base = NAN, htop = NAN;
    float meanu = NAN, meanv = NAN, uchg = 0, vchg = 0, srmag = 0;

    int i=0,k=0;

    if ( (*mu).cape > 100.0 ){/*parcel method*/

        /*estimate bottom of the effective inflow layer*/
        while (i < anzahl){
            if ( (all[i].cape >= 100) && (all[i].cin > -250) ){
                //pbot =  (float) drk[n].druck;
                base =  (float) all[i].agl;
                break;
            }
            i++;
        }

        depth = (*mu).eqh - base;
        htop = base + ( depth * 0.65 );

        meanu = mean_wind_u(all,base,htop);
        meanv = mean_wind_v(all,base,htop);
        k=i;
        while ( (k < anzahl) && (all[k].agl < htop) ){
            k++;
        }
        printf("Parcel Method: %.1f\t%.1f\n",htop,all[k].agl);
        shru = all[k].u - all[i].u;
        shrv = all[k].v - all[i].v;   

        srmag = sqrt( SQR(shru) + SQR(shrv) );
        uchg = d / srmag * shrv;
        vchg = d / srmag * shru;

        (*bu).rstu = meanu + uchg;
        (*bu).rstv = meanv - vchg;
        (*bu).lstu = meanu - uchg;
        (*bu).lstv = meanv + vchg;

    } else {/*non parcel method*/
        // SFC-6km Shear Vector
        // Experimentel
        shru = mean_wind_u(all,5500,6000) - mean_wind_u(all,0,500);
        shrv = mean_wind_v(all,5500,6000) - mean_wind_v(all,0,500);
        // normal method
        // don't forget to uncomment the variables
        /*
        while (drk[i].agl < 6000){
            i++;
        }  
        //wmean
        dz = drk[i].agl - drk[i-1].agl;
        u_top = drk[i].u * ((high - drk[i-1].agl)/dz) +  drk[i-1].u * ((drk[i].agl - high)/dz);
        v_top = drk[i].v * ((high - drk[i-1].agl)/dz) +  drk[i-1].v * ((drk[i].agl - high)/dz);

        shru = u_top - drk[0].u;
        shrv = v_top - drk[0].v;
        */
        meanu6 = mean_wind_u(all,0,6000);
        meanv6 = mean_wind_v(all,0,6000);

        tmp = d / sqrt( SQR(shru) + SQR(shrv) );
        
        (*bu).rstu = meanu6 + (tmp * shru);
        (*bu).rstv = meanv6 - (tmp * shrv);
        (*bu).lstu = meanu6 - (tmp * shru);
        (*bu).lstv = meanv6 + (tmp * shrv);
//        printf("\n %.1f\n%.1f\n",meanu6,meanv6);
//        printf("\n%.1f\n%.1f\n",shru,shrv);

    }
    printf("Storm Motion: rstu:%.1f,rstv:%.1f,lstu:%.1f,lstv%.1f\n\n",(*bu).rstu,(*bu).rstv,(*bu).lstu,(*bu).lstv);
    return bu;
    free (bu);
}
/*End Storm Motion Vector*/


/*Storm relative winds*/

void stormrelativwind (struct stormmotion *ptr1, struct sounding_data *ptr2){
    /*Wind relative to a spuercell, which follow the storm motion vector
     *Diffrence between u,v components and storm motion vectors*/

    (*ptr2).srur = (*ptr2).u - (*ptr1).rstu;
    (*ptr2).srvr = (*ptr2).v - (*ptr1).rstv;
    (*ptr2).srul = (*ptr2).u - (*ptr1).lstu;
    (*ptr2).srvl = (*ptr2).v - (*ptr1).lstv;
    return;
}


float srh ( struct sounding_data drk[], float high, int mode){

    float srh=0,h=0;
    int i=0;

    while (drk[i].agl < high){
        if (mode == 1){
            //srh += ( ( drk[i].ushr * drk[i].srvr ) - (drk[i].vshr * drk[i].srur ));
            srh +=  drk[i].ushr * ((drk[i].srvr + drk[i+1].srvr)/ (float) 2) - drk[i].vshr * ((drk[i].srur +  drk[i+1].srur)/ (float) 2);
        } else if (mode == 2){
            //srh += ( drk[i].ushr * drk[i].srvl ) - (drk[i].vshr * drk[i].srul );
            srh +=  drk[i].ushr * ((drk[i].srvl + drk[i+1].srvl)/ (float) 2) - drk[i].vshr * ((drk[i].srul +  drk[i+1].srul)/ (float) 2);
        } else if (mode == 3){
            if( ( h = drk[i].ushr * ((drk[i].srvr + drk[i+1].srvr)/ (float) 2) - drk[i].vshr * ((drk[i].srur +  drk[i+1].srur)/ (float) 2) ) >= 0)
                srh += h;
        } else {
            srh += 0;
        }
        i++;
    }

//    if (mode == 2)
//        srh *= -1;

    return srh;
}


/*Effective Inflow Layer*/

/*
struct inflow {
    float ptop;
    float pbot;
    float esrh;
    float ebwd;
};
*/

//effective_inflow_layer(prof, 100, -250, mupcl=mupcl)
//temp.efi = (*effective_inflow_layer(drk,temp.mu.eqh,j));


struct inflow *effective_inflow_layer (struct sounding_data drk[],int eqh,int anzahl) {

    struct inflow *efi;
    efi = malloc (sizeof(struct inflow));
    int i=0,index=-1; 
    (*efi).esrh=0.0;

    /*Find bottom layer for potential thunderstorm updraft*/
    while ( i < anzahl ){
        if ( ( drk[i].cape >= 100) && ( drk[i].cin >= -250 ) ){
            (*efi).pbot = drk[i].druck;
            index = i;
            break;
        }
        i++;
    }

    /*Find top layer for potential thunderstorm updraft*/
    while (i < anzahl){
        if ( ( drk[i].cape >= 100) && ( drk[i].cin >= -250 ) ){
            (*efi).esrh += drk[i].ushr * ((drk[i].srvr + drk[i+1].srvr)/(float) 2) - drk[i].vshr * ((drk[i].srur + drk[i+1].srur)/ (float) 2); 
            i++;
        } else {//hier auch wmean
            i--;
            (*efi).ptop = drk[i].druck;
            break;
        }
    }

    i = index;
    eqh /= 2;
    while ( (int) drk[i].agl < eqh) {
        i++;
    }
    if (index < 0 ) {
        (*efi).esrh = 0;
        (*efi).ebwd = 0;
        (*efi).pbot = NAN;
        (*efi).ptop = NAN;
    } else {
        /*Calculate ebwd*/
        //(*efi).esrh = esrh;
        //float shear (float ubot, float vbot, float utop, float vtop)
        (*efi).ebwd = shear (drk[index].u,drk[index].v,drk[i].u,drk[i].v);
    }
    printf("\nMessage from Effective Inflow Layer:\n\tpbot : %.1f\n\tptop : %.1f\n\tebwd : %.1f\n\tesrh : %.1f\n",(*efi).pbot,(*efi).ptop,(*efi).ebwd ,(*efi).esrh);
    return efi;
    free(efi);
}


float streamwise (struct sounding_data all[], float high) {

    int i=0;
    float H=0, h=0, st=0;

    while (all[i].agl < high){
        h =  all[i].ushr * ((all[i].srvr + all[i+1].srvr)/ (float) 2) - all[i].vshr * ((all[i].srur +  all[i+1].srur)/ (float) 2);
        st = sqrt( SQR((all[i].srur +  all[i+1].srur)/ (float) 2) + SQR((all[i].srvr + all[i+1].srvr)/ (float) 2));
//      printf("%.1f\n",h/st );
        if (st != 0)
            H += h/st;
        i += 1;
        //printf("%.1f",h);
    }
    return H;
}


float critical_angle (struct indices *ptr1, struct sounding_data all[], float high) {

    int i=0;
    float u500=0, v500=0, sfc_u=0, sfc_v=0,dz=0,dot=0,angle=0;
    float vec1_u=0, vec1_v=0, vec2_u=0, vec2_v=0, vec_1_mag=0, vec_2_mag=0;
    float val = 180 / PI;

    sfc_u=all[0].u;
    sfc_v=all[0].v;

    while (all[i].high < high){
        i++;
    }
    //wmean
    dz = all[i].high - all[i-1].high;
    u500 = all[i].u * ((high - all[i-1].high)/dz) + all[i-1].u * ((all[i].high - high)/dz);
    v500 = all[i].v * ((high - all[i-1].high)/dz) + all[i-1].v * ((all[i].high - high)/dz);
    printf("SFC Wind: %.1f\t%.1f\n",sfc_u,sfc_v);

    vec1_u = u500 - sfc_u;
    vec1_v = v500 - sfc_v;
    vec2_u = (*ptr1).bunker.rstu - sfc_u;
    vec2_v = (*ptr1).bunker.rstv - sfc_v;
//    vec2_u = sfc_u;
//    vec2_v = sfc_v;
    vec_1_mag = sqrt(SQR(vec1_u) + SQR(vec1_v) );
    vec_2_mag = sqrt(SQR(vec2_u) + SQR(vec2_v) );

    dot = vec1_u * vec2_u + vec1_v * vec2_v;
    angle = val*(acos(dot / (vec_1_mag * vec_2_mag)));
    return angle;
}


float length(struct sounding_data all[], float bot, float top){
    /*
    Calculate the length of the hodograph.

    Paramter:
    ---------
    ushr : Shear
    vshr : Shear
    height : till height
    bot : in m
    top : in m

    Returns:
    --------
    length in m/s
    
    
    */
    float length = 0;
    int i = 0;
    while( (all[i].agl >= bot) && (all[i].agl < top) ){
        length += sqrt(SQR( all[i].ushr ) + SQR( all[i].vshr ));
        i++;
    }

    return length;
}
/*
float lemgth_major_layers(struct sounding_data all[]){

    float length = 0;
    float height[] = [1000,925,850,700,500,400,300]
    while(all[i].druck == height){
         
    }
}
*/
#endif /*SHEAR_H*/
