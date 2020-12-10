#include<math.h>
#include<stdio.h>
#include<stdlib.h>

/* Random number generator */

#define PI 3.14159265358979323846
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* End of random number definition */

/* In the model, the lengthscale is in nm, the timescale is in second, the force is in pN */

/* declare of global variables */

double K_A=25.6*0.5;            /* the spring constant of local binding potential by the FtsZ monomer, the unit is pN/nm, 25.6 corresponds to the potential well of 20 KbT with the range of +/-2.5 nm*/
double D_I=40000.0;             /* the intrinisic diffusion constant of FtsI, the unit is nm^2/second */
double tau_Z=0.2;               /* the lifetime of FtsZ monomer at the end of the FtsZ filament, in the unit of second */
double tau_N=0.6;               /* the lifetime of FtsN on the denuded PG strand, in the unit of second */
double dt=0.000005;             /* simulation time step, in the unit of second */
int N_Z=40000;                  /* the number of time steps that FtsZ monomer will fall off, 20000 steps correspond to 0.2 second */
int N_A=40000;                  /* the number of time steps that FtsA will fall off, together with the FtsZ monomer below, 100000 steps correspond to  1 second */
int N_N=2000000;                /* the number of time steps that FtsN on the leftmost will fall off, 60000 steps correspond to 0.6 second */
int T_max=12000000;             /* the maximum number of simulation time steps */
int T_interval=11999999;        /* the interval of simulation time steps for output */
int Z_length=51;                /* FtsZ filament length – number of subunits */

/* end of declaration */

/* the seed for random number generator */

long *idumA;
long *idumB;

/* declare of functions */

float ran1A(long *idumA);      
float gasdevA(long *idumA);      /* random number generator for the FtsI position */

float ran1B(long *idumB);
float gasdevB(long *idumB);

/* end of declare of functions */


main(){

    FILE *positions_vs_time;

    double noise_level_A=sqrt(2.0*D_I*dt);       /* noise level of the FtsI position at each simulation time step, in the unit of nm.  It equals to the effect of diffusion.  */
    double zeta=4.0/D_I;                         /* the viscous drag coefficient of FtsI, 1 KbT/D_I, in the unit of pN.second/nm */
    
    double X_I=0.0;                /* The initial position of FtsI, 255 nm = 51* 5nm */
    double F_Z=0.0;                /* The force on FtsI by FtsZ */
    double F_A=0.0;                /* The force on FtsI by FtsA */
    
    int t=0;                           /* The index for simulation steps */
    double X_Z[1200], X_A[1200];      /* The positions of FtsZ and FtsA monomers */
    double X_Z1[1200], X_A1[1200];
    int State_Z[1200], State_A[1200]; /* The on and off states for FtsZ and FtsA monomers */
    int i, j, k;                      /* The indexes for FtsZ and FtsA monomers */
    int end_Z=1;                   /* The index for FtsZ shrinking end, initially is 1, will increase 1 by every tau_Z second, which is 0.2 second and 2000 simulation steps */
    int end_A=1;                   /* The index for FtsA/Z shrinking end, initially is 1, will increase 1 by every 5*tau_Z second, which is 1 second and 10000 simulation steps */

    int end_AR=Z_length;           /* The index for FtsA/Z growing end, initially is 51*/
    
    double lower_boundary=0.0;
    
    double tracking_TL=0.0;        /* recording the time that FtsI is < 100nm away from FtsZ leftmost end */
    double tracking_TR=0.0;        /* recording the time that FtsI is < 100nm away from FtsZ rightmost end */
    int tracking_T_first_encounter=0;   /* recording the first time that FtsI is < 100nm within the FtsZ filament tips */
    int record_first_encouter=0;
    double tracking_T_first_encounter_position=0.0;     /* recording the position that FtsI is <100 nm within the FtsZ filament tips for the first time */
    double total_time_bound=0.0;        /* record the total number of steps that FtsI is bound to FtsZ, +/- 100nm from both ends */
    
    
    
    int temp_Z=0;
    int temp_A=0;

    int temp_Z1=0;               /* temporary position of the other FtsI along Z/A filament */
    int temp_A1=0;
    
    
    int temp1=0;
    int temp2=0;

    int mark=0;

    int timer=0;                    /* internal clock recording the lifetime of active FtsI */
    int temp_timer=0;               /* temporary index equals to the instantaneous "timer" */

    long dummyA=-5;
    long dummyB=-100;                /* random number initiation for FtsZ treadmilling speed */

    idumA=&dummyA;
    idumB=&dummyB;
    

    for(i=0;i<=1199;i++){    /* Initialization of FtsZ monomer positions and the corresponding state. State_Z[i]=1 means "on", and 0 means "off" */
        X_Z[i]=5.0*i;
        X_Z1[i]=X_Z[i]-100.0;
        if(i==0){
            State_Z[i]=0;
        }
        else{
            State_Z[i]=1;
        }
    }

    for(j=0;j<=1199;j++){    /* Initialization of FtsA monomer positions and the corresponding state. State_A[i]=1 means "on", and 0 means "off" */
        X_A[j]=5.0*j;
        X_A1[j]=X_A[j]-100.0;
        if(j==0){
            State_Z[j]=0;
        }
        else{
            State_Z[j]=1;
        }
    }

    
    for(k=0;k<1;k++){
     
        X_I=5.0*ran1B(idumB);                  /* FtsI initial position is random between 0 and 5 nm;  if it times "Z_length" will render the FtsI initial position between 0 and 255 nm */

        int N_Z_ran=N_Z*(1+0.3*gasdevB(idumB));    /* the number of time steps it takes for Z/A end monomer to fall off in subject to Gaussian fluctutations, with 30% of SD */
        int N_A_ran=N_A*(1+0.3*gasdevB(idumB));

        
        tracking_TL=0.0;
        tracking_TR=0.0;
        tracking_T_first_encounter=0;
        tracking_T_first_encounter_position=0.0;
        record_first_encouter=0;
        total_time_bound=0.0;
        
        end_A=1;
        end_Z=1;
        end_AR=Z_length;
        
        F_A=0.0;
        F_Z=0.0;

        for(t=0;t<T_max;t++){
            
            if(fabs(X_I-X_A[end_A])<100.0){     /* record the time that FtsI end-tracks the shrinking end of FtsZ filament */
                tracking_TL+=1.0;
            }
            
            if(fabs(X_I-X_A[end_AR])<100.0){    /* record the time that FtsI end-tracks the growing end of FtsZ filament */
                tracking_TR+=1.0;
            }
            
            if(X_I>=(X_A[end_A]-100.0)&&X_I<=(X_A[end_AR]+100.0)){    /* record the time that FtsI within the length of FtsZ filament +/- 100 nm */
                total_time_bound+=1.0;
            }
            
            if(((X_A[end_A]-X_I)>100.0||(X_I-X_A[end_AR])>100.0)&&(record_first_encouter==0)){   /* record the first time that the FtsI drops from the FtsZ shrinking end more than 100 nm */
                tracking_T_first_encounter=t;
                tracking_T_first_encounter_position=X_I;
                record_first_encouter=1;
            }
            
            /* if the FtsI does not fall off the track throughout the 60-sec, then its first-passage-time of escape is set to be the length of simulation time */
            
            if((t==(T_max-1))&&(tracking_T_first_encounter==0)){
                tracking_T_first_encounter=T_max-1;
                tracking_T_first_encounter_position=X_I;
            }
            
            
            if(t%T_interval==0&&t>0){
                positions_vs_time=fopen("/Users/liuj7/CPrograms/FtsZ_treadmilling/output/XZ_XA_XI_XN_time.txt", "a");
            }
            
            end_Z=(int)(t/N_Z_ran)+1;    /* The index for the leftmost FtsZ monomer that still exists with the "on" state */
            end_A=(int)(t/N_A_ran)+1;    /* The index for the leftmost FtsZ monomer that still exists with the "on" state */
            
            end_AR=(int)(t/N_A_ran)+Z_length;   /* The index for the rightmost FtsZ monomer that still exists with the "on" state */
            
            temp_Z=(int)(X_I/5.0);                /* determine the index of the corresponding FtsZ monomer that is directly underneath the FtsI */
            temp_A=(int)(X_I/5.0);                /* determine the index of the corresponding FtsA monomer that is directly underneath the FtsI */
            
            
            if(end_Z>temp_Z){                     /* if the FtsI is outside the interior of FtsZ filament */
                if((X_Z[end_Z]-X_I)>2.5){
                    F_Z=0.0;
                }
                else{
                    F_Z=X_Z[end_Z]-X_I;
                }
                
            }
            else{
                F_Z=-X_I+5.0*temp_Z;
            }
            
            
            if(end_A>temp_A){
                if((X_A[end_A]-X_I)>2.5){
                    F_A=0.0;
                }
                else{
                    F_A=X_A[end_A]-X_I;
                }
            }
            
            else if(end_AR<temp_A){
                if((X_I-X_A[end_AR])>2.5){
                    F_A=0.0;
                }
                else{
                    F_A=X_A[end_AR]-X_I;
                }
            }
            
            else{
                if((X_I-X_A[temp_A])>2.5){
                    F_A=5.0*temp_A-X_I+5.0;
                }
                else{
                    F_A=5.0*temp_A-X_I;
                }
            }
            
            X_I+=(K_Z*F_Z+K_A*F_A)*dt/zeta+noise_level_A*gasdevA(idumA);
            
            if(t%T_interval==0&&t>0){
                fprintf(positions_vs_time, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", 1.0*t*dt, X_I, X_A[end_A], X_A[end_AR], 5.0/(N_A_ran*dt), tracking_TL/T_max, tracking_TR/T_max, total_time_bound*dt, tracking_T_first_encounter*dt, tracking_T_first_encounter_position/(tracking_T_first_encounter*dt+0.00001));
            }
            
            if(t%T_interval==0&&t>0){
                fclose(positions_vs_time);
            }
        }
    }
}



/* Random number generator 0-1 */
float ran1A(long *idumA)                 
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idumA<=0||!iy) {
    if(-(*idumA)<1) *idumA=1;
    else *idumA=-(*idumA);
    for(j=NTAB+7;j>=0;j--){
      k=(*idumA)/IQ;
      *idumA=IA*(*idumA-k*IQ)-IR*k;
      if(*idumA<0) *idumA+=IM;
      if(j<NTAB) iv[j]=*idumA;
    }
    iy=iv[0];
  }
  k=(*idumA)/IQ;
  *idumA=IA*(*idumA-k*IQ)-IR*k;
  if(*idumA<0) *idumA+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=*idumA;
  if((temp=AM*iy)>RNMX) return RNMX;
  else return temp;
}


/* Gaussian distribution, mean 0, stdev 1, using ran1(idum) */

float gasdevA(long *idumA)
{
  float ran1A(long *idumA);
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;

  if(iset==0){
    do{
      v1 = 2.0*ran1A(idumA)-1.0;
      v2 = 2.0*ran1A(idumA)-1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq>=1.0||rsq==0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);

    gset=v1*fac;
    iset= 1;
    return v2*fac;
  }
  else{
    iset = 0;
    return gset;
  }
}


/* Random number generator 0-1 */
float ran1B(long *idumB)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    if (*idumB<=0||!iy) {
        if(-(*idumB)<1) *idumB=1;
        else *idumB=-(*idumB);
        for(j=NTAB+7;j>=0;j--){
            k=(*idumB)/IQ;
            *idumB=IA*(*idumB-k*IQ)-IR*k;
            if(*idumB<0) *idumB+=IM;
            if(j<NTAB) iv[j]=*idumB;
        }
        iy=iv[0];
    }
    k=(*idumB)/IQ;
    *idumB=IA*(*idumB-k*IQ)-IR*k;
    if(*idumB<0) *idumB+=IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j]=*idumB;
    if((temp=AM*iy)>RNMX) return RNMX;
    else return temp;
}


/* Gaussian distribution, mean 0, stdev 1, using ran1(idum) */

float gasdevB(long *idumB)
{
    float ran1B(long *idumB);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    
    if(iset==0){
        do{
            v1 = 2.0*ran1B(idumB)-1.0;
            v2 = 2.0*ran1B(idumB)-1.0;
            rsq = v1*v1 + v2*v2;
        } while (rsq>=1.0||rsq==0.0);
        fac = sqrt(-2.0*log(rsq)/rsq);
        
        gset=v1*fac;
        iset= 1;
        return v2*fac;
    }
    else{
        iset = 0;
        return gset;
    }
}
