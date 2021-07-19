
//#1 DECLARATION SEQUENCE------------------------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
//#include <conio.h>
#include <stdlib.h>
#include <stdio.h>


//R_A_N_2_P_A_R_A_M_E_T_E_R_S___________________________________________________

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 3.14159


//_P_A_R_A_M_E_T_E_R_S__________________________________________________________

#define T 2500//Time of simulation in msec

#define S 1//size of parameter space
#define Trials 10//number of trials (if averaging is needed)

#define Ne 800//880//480 //number of excitatory units
#define Ni 200//220//120// number of  inhibitory units
#define dt 0.1 //time step for integration
#define Dt 0.001// Value in sec of 1 dt e.g. Dt=0.001 means that 1 dt equals 1 msec

#define temporal_window 20 //Temporal window for spike coherence measures

#define beta 300 //300//activation function non-linear gain

#define gamma 0.016 //scaling to mV

#define threshold  0.0//activation function threshold (inflexion point)

#define p_connectivity 1.00//0.2 //probability of synaptic connexion

#define D0 0.001// Noise variance

#define I_e0 -0.25//Bias cuurent e cells !ADD CURRENT HERE TO TRIGGER SEIZURE LIKE STATE/GAMMA
//#define I_i0 -0.5// Bias current i cells
#define I_i0 -.5// Bias current i cells
#define c 0.0 //noise spatial correlations

//#define sigma_e 0.05//heterogeneity in e cells VARIED BELOW THATS WHY THESE ARE COMMENTED OUT HERE
//#define sigma_i 0.05//heterogeneity in e cells

#define w_ee0  1.6   //e->e synaptic strenght
#define w_ie0  -4.7 //i->e synaptic strenght
#define w_ei0  3.0   //e->i synaptic strenght
#define w_ii0  -.13 //i->i synaptic strenght

#define alpha_e 1 //time scale E units
#define alpha_i 2 //time scale I units

#define b 0// adaptation gain
#define alpha_adapt 0.01//adaptation time scale - should be slower

double f(double u, double h); // Firing rate function
double F(double u, double sigma);
double f_prime(double u,double h);
double f_2prime(double u,double h);
float ran2(long *idum);//random number generator - initial cond.
void shuffle_indices(int Array[], int size);//shuffle an array with int  entries
void shuffle(double Array[], int size);//shuffle an array with double  entries
void sort(double Array[], int size);//sort array 
void four1(double data[], unsigned long nn, int isign);
double coherence(double X1[], double X2[], int window);

double PSD[(int)T];
double Freq[(int) T];
double Mean_PSD[(int)T];
int order_e[Ne];
int order_i[Ni];
double u_e[Ne][T];//Dynamics of E units
double u_i[Ni][T];//Dynamics of I units
double U_e[T];
double V_e[T];
double U_i[T];
double V_i[T];
double a_e[Ne][T];
double a_i[Ni][T];
double mean_u_e[T];//mean network activity E units
double mean_u_i[T];//mean network activity I units
double Input[T];
double rate_e[T];
double rate_i[T];
double X_e[Ne][T];
double X_i[Ni][T];//spikes
double Xi_e[Ne][T];
double Xi_i[Ni][T];
double Xi_c[T];
double h_e[Ne];
double h_i[Ni];
double eta_e[Ne];
double eta_i[Ni];
double Spikes_e[Ne][T];
double Spikes_i[Ne][T];
double COHERENCE[S][S];
double Sigma_e[S];
double Sigma_i[S];
double Indiv_Rates_e[Ne];
int Display_indices_e[Ne];
using namespace std;

double connectivity_ee[Ne][Ne];
double connectivity_ii[Ni][Ni];
double connectivity_ei[Ni][Ne];
double connectivity_ie[Ne][Ni];

char output1[100];
char output2[100];
char output3[100];
char output4[100];
char output5[100];


int main (int argc, char *argv[])
{
	double sigmae;
	double sigmai;

	sigmae=atof(argv[1]);
	sigmai=atof(argv[2]);
    
    srand (time(NULL));
    
    int n;
    ofstream outfile;
    n=sprintf(output3,"SpikingNEWNEWv2_p%2.0f_sige%2.0f_sigi%2.0f_Trials%d_gamma%3.0f.csv", p_connectivity*100, sigmae*10, sigmai*10, Trials, gamma*1000);
    outfile.open(output3, ios::out);
    
    ofstream outfile3;
    n=sprintf(output4,"SpikingInhibNEWNEWv2_p%2.0f_sige%2.0f_sigi%2.0f_Trials%d_gamma%3.0f.csv", p_connectivity*100, sigmae*10, sigmai*10, Trials, gamma*1000);
    outfile3.open(output4, ios::out);
    
    ofstream outfile2;
    n=sprintf(output2,"DynamicsNEWNEWv2_p%2.0f_sige%2.0f_sigi%2.0f_Trials%d_gamma%3.0f.csv", p_connectivity*100, sigmae*10, sigmai*10, Trials, gamma*1000);
    outfile2.open(output2, ios::out);

    ofstream outfile4;
    n=sprintf(output4,"hNEWNEWv2_p%2.0f_sige%2.0f_sigi%2.0f_Trials%d_gamma%3.0f.csv", p_connectivity*100, sigmae*10, sigmai*10, Trials, gamma*1000);
    outfile4.open(output4, ios::out);

    ofstream outfile5;
    n=sprintf(output5,"hInhibNEWNEWv2_p%2.0f_sige%2.0f_sigi%2.0f_Trials%d_gamma%3.0f.csv", p_connectivity*100, sigmae*10, sigmai*10, Trials, gamma*1000);
    outfile5.open(output5, ios::out);
    
    double w_ee;
    double w_ei;
    double w_ie;
    double w_ii;
    double D;
    double I_i;
    double I_e;
    
    w_ee=w_ee0/gamma;
    w_ei=w_ei0/gamma;
    w_ie=w_ie0/gamma;
    w_ii=w_ii0/gamma;
    D=D0/(gamma*gamma);
    I_i=I_i0/gamma;
    I_e=I_e0/gamma;

    cout<<w_ee<<"    "<<w_ei<<"    "<<w_ie<<"    "<<w_ii<<"    "<<D<<"    "<<I_i<<endl;

    
    
    for (int s1=0;s1<.5;s1++)
    {
        
        
        for (int s2=0;s2<.5;s2++)
        {
            
            
//            Sigma_e[s1]=0.05+0.25*s1/((double)S);
//            Sigma_i[s2]=0.05+0.25*s2/((double)S);
            
            Sigma_e[s1]=sigmae;
            Sigma_i[s2]=sigmai;

            
            
            //initial conditions
            for (int t=0;t<T;t++)
            {
                Mean_PSD[t] =0;
                
                for(int n=0;n<Ne;n++)
                {
                    u_e[n][t]=0;
                    rate_e[n]=0;
                    Xi_e[n][t]= 0;
                    X_e[n][t] = 0;
                    order_e[n]=0;
                    
                }
                for(int n=0;n<Ni;n++)
                {
                    u_i[n][t]=0;
                    rate_i[n]=0;
                    Xi_i[n][t]= 0;
                    X_i[n][t] = 0;
                    order_i[n]=0;
                }
                mean_u_e[t]=0;
                mean_u_i[t] =0;
                mean_u_e[t]=0;
                PSD[t]=0;
                Xi_c[t]=0;
                
            }
            
            
            
            for (int q=0;q<Trials;q++)
            {
                
                
                cout<<s1<<"	"<<s2<<"	"<<q<<endl;
                //clean memory
                for (int t=0;t<T;t++)
                {
                    for(int n=0;n<Ne;n++)
                    {
                        u_e[n][t]=0;
                        Xi_e[n][t]= 0;
                        X_e[n][t] = 0;
                        order_e[n]=0;
                        a_e[n][t]=0;
                        
                        
                    }
                    for(int n=0;n<Ni;n++)
                    {
                        u_i[n][t]=0;
                        Xi_i[n][t]= 0;
                        X_i[n][t] = 0;
                        order_i[n]=0;
                        
                    }
                    mean_u_e[t]=0;
                    mean_u_i[t] =0;
                    mean_u_i[t]=0;
                    PSD[t]=0;
                    Input[t]=0;
                    
                    for(int n=0; n<Ne; n++)
                    {
                        for(int m=0; m<Ne;m++)
                        {
                            connectivity_ee[n][m]=0;
                        }
                    }
                    
                    for(int n=0; n<Ni; n++)
                    {
                        for(int m=0; m<Ne;m++)
                        {
                            connectivity_ei[n][m]=0;
                        }
                    }
                    
                    for(int n=0; n<Ni; n++)
                    {
                        for(int m=0; m<Ni;m++)
                        {
                            connectivity_ii[n][m]=0;
                        }
                    }
                    
                    for(int n=0; n<Ne; n++)
                    {
                        for(int m=0; m<Ni;m++)
                        {
                            connectivity_ie[n][m]=0;
                        }
                    }
                }
                
                
                // Create connectivity matricies
                for(int n=0; n<Ne; n++)
                {
                    for(int m=0; m<Ne; m++)
                    {
                        if (m!=n)
                        {
                            if (((double) rand()/ RAND_MAX)<p_connectivity)
                            {
                                connectivity_ee[n][m]=1;
                            }
                            
                        }
                    }
                }
                
                
                for(int n=0; n<Ni; n++)
                {
                    for(int m=0; m<Ne; m++)
                    {
                        if (m!=n)
                        {
                            if (((double) rand()/ RAND_MAX)<p_connectivity)
                            {
                                connectivity_ei[n][m]=1;
                            }
                            
                        }
                    }
                }
                
                
                for(int n=0; n<Ni; n++)
                {
                    for(int m=0; m<Ni; m++)
                    {
                        if (m!=n)
                        {
                            if (((double) rand()/ RAND_MAX)<p_connectivity)
                            {
                                connectivity_ii[n][m]=1;
                            }
                            
                        }
                    }
                }
                
                
                for(int n=0; n<Ne; n++)
                {
                    for(int m=0; m<Ni; m++)
                    {
                        if (m!=n)
                        {
                            if (((double) rand()/ RAND_MAX)<p_connectivity)
                            {
                                connectivity_ie[n][m]=1;
                            }
                            
                        }
                    }
                }
                //                for(int n=0; n<Ni; n++)
                //                {
                //                cout<<connectivity_ii[0][n];
                //                }
                //                cout<<" "<<endl;
                
                
                
                
                
                
                
                
                
                // define noise arrays AND INPUT
                for (int t=0;t<T;t++)
                {
                    long d=rand()%105;
                    long noisyseed5=(long)15*t+3+8+q*44+s2*d*78+34*s2;
                    long noisyseed6=(long)9*t+t*3+q*1+56*d+341*s2+1;
                    Xi_c[t]= sqrt(-2*log(ran2(&noisyseed5)))*cos(2*3.14159*ran2(&noisyseed6));
                    for(int n=0;n<Ne;n++)
                    {
                        long d=rand()%105;
                        long noisyseed1=(long)21*t+n+q*187198+56*d+12*s1+1;
                        long noisyseed2=(long)69*t*q+11*n+t+45*d+s2*s1+s1;
                        
                        Xi_e[n][t]= sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
                        X_e[n][t] = 0;
                        
                        
                        
                        
                    }
                    for(int n=0;n<Ni;n++)
                    {
                        long d=rand()%87;
                        long noisyseed3=(long)15*t+3+n*8+q*8976+s2*d*78+34*s2;
                        long noisyseed4=(long)9*t+3*n+t*3+q*38493+56*d+341*s2+1;
                        Xi_i[n][t]= sqrt(-2*log(ran2(&noisyseed3)))*cos(2*3.14159*ran2(&noisyseed4));
                        X_i[n][t] = 0;
                        
                    }
                    
                    {
                        Input[t]=(double)((double)((double)(t)/(double)(T))*.5)/gamma;
//                        cout<<Input[t]<<q<<endl;

                    }
                    
//                    if (t>1095 && t<1105)
//                    {
//                        Input[t]=-I_e*2;
//                    }
//
//                    if (t>1595 && t<1605)
//                    {
//                        Input[t]=-I_e*2;
//                    }
                    
                }
                
                //sample rheobases
                for (int k=0;k<Ne;k++)
                {
                    long d=rand()%55;
                    long noisyseed7=(long)56*d+k+s1*s2+67*s2;
                    long noisyseed8=(long)556*d+15*k+s1*s2+671*s2;
                    eta_e[k]= sqrt(-2*log(ran2(&noisyseed7)))*cos(2*3.14159*ran2(&noisyseed8));
                    h_e[k] = (Sigma_e[s1])*eta_e[k];
        			outfile4<<h_e[k]<<endl;             
                }
                for (int k=0;k<Ni;k++)
                {
                    long d=rand()%85;
                    long noisyseed9=(long)88*d+k*89+12*s1*s2+6*s2*88;
                    long noisyseed10=(long)666*d+12*k+121*s1*s2+333*s2;;
                    eta_i[k]= sqrt(-2*log(ran2(&noisyseed9)))*cos(2*3.14159*ran2(&noisyseed10));
                    h_i[k] = (Sigma_i[s2])*eta_i[k];
        			outfile5<<h_i[k]<<endl;
                }
                
                
                
                
                
                
                //Initial conditions
                for (int t=0;t<T;t++)
                {
                    for (int n=0;n<Ne;n++)
                    {
                        u_e[n][t]=Xi_e[n][t];
                        
                        
                    }
                    for (int n=0;n<Ni;n++)
                    {
                        u_i[n][t]=Xi_i[n][t];
                        
                    }
                    
                }
                
                
                
                
                
                
                
                
                ////////////////////EULER SCHEME//////////////////////////////////////////////////////////////////////////
                
                //Define integration order
                
                for (int i=0;i<Ne;i++)
                {
                    order_e[i]=i;
                }
                for (int i=0;i<Ni;i++)
                {
                    order_i[i]=i;
                }
                shuffle_indices(order_e,Ne);
                shuffle_indices(order_i,Ni);
                
                
                double mean_ue=0;
                double mean_ui=0;
                double mean_p_e=0;
                double mean_p_i=0;
                double mean_ee=0;
                double mean_ei=0;
                double mean_ie=0;
                double mean_ii=0;
                
                // Integrate network dynamics
                
                for (int t=0;t<T;t++)
                {
                    shuffle_indices(order_e,Ne);
                    shuffle_indices(order_i,Ni);
                    
                    
                    
                    
                    for (int i=0;i<Ne;i++)
                    {
                        
                        double sum_ee=0;
                        
                        
                        
                        
                        for (int j=0;j<Ne;j++)
                        {
                            
                            sum_ee=sum_ee+1/((double)Ne)*(w_ee/p_connectivity)*connectivity_ee[i][j]*X_e[j][t];
                        }
                        mean_ee=mean_ee+1/((double)T)*1/((double)Ne)*sum_ee;
                        
                        
                        
                        double sum_ie=0;
                        
                        
                        for (int j=0;j<Ni;j++)
                        {
                            sum_ie=sum_ie+1/((double)Ni)*(w_ie/p_connectivity)*connectivity_ie[i][j]*X_i[j][t];
                            
                            
                        }
                        mean_ie=mean_ie+1/((double)T)*1/((double)Ni)*sum_ie;
                        
                        
                        
                        u_e[order_e[i]][t+1]=(u_e[order_e[i]][t]+dt*alpha_e*(-1*u_e[order_e[i]][t]+b*a_e[order_e[i]][t]+sum_ee+sum_ie+I_e+1*Input[t+1])+sqrt(2*alpha_e*D*dt)*(sqrt(1-c)*Xi_e[order_e[i]][t]+sqrt(c)*Xi_c[t]));
                        a_e[order_e[i]][t+1]= a_e[order_e[i]][t]+dt*alpha_adapt*(1*(-a_e[order_e[i]][t]+(u_e[order_e[i]][t]-I_e)));
                        
                        //Spiking activity of excitatory neurons
                        
                        // Poisson neurons
                        
                        long seed_q =rand()%101;
                        long seed_e = (long) 12+i*i+t+15+t*seed_q+56*s1+s2+99*q+2;
                        double p_e =ran2(&seed_e);
                        double p_fire_e = ( 1-exp(-f(u_e[order_e[i]][t+1],h_e[order_e[i]])*dt));
                        
                        
                        if (p_e<p_fire_e)//spike occurs
                            
                        {
                            X_e[order_e[i]][t+1] =1/dt;
                            
                            
                        }
                        else//no spike
                        {
                            X_e[order_e[i]][t+1]=0;
                            
                        }
                        
                        
                        
                    }
                    
                    
                    
                    for (int i=0;i<Ni;i++)
                    {
                        double sum_ei=0;
                        
                        
                        for (int j=0;j<Ne;j++)
                        {
                            
                            sum_ei=sum_ei+1/((double)Ne)*(w_ei/p_connectivity)*connectivity_ei[i][j]*X_e[j][t];
                            
                            
                        }
                        
                        double sum_ii=0;
                        
                        
                        for (int j=0;j<Ni;j++)
                        {
                            
                            sum_ii=sum_ii+1/((double)Ni)*(w_ii/p_connectivity)*connectivity_ii[i][j]*X_i[j][t];
                        }
                        mean_ii=mean_ii+1/((double)T)*1/((double)Ni)*sum_ii;
                        
                        
                        u_i[order_i[i]][t+1]=(u_i[order_i[i]][t]+dt*alpha_i*(-u_i[order_i[i]][t]+b*a_i[order_i[i]][t]+sum_ei+sum_ii+I_i+0*Input[t+1])+sqrt(2*alpha_i*D*dt)*(sqrt(1-c)*Xi_i[order_i[i]][t]+sqrt(c)*Xi_c[t]));
                        a_i[order_i[i]][t+1]= a_i[order_i[i]][t]+dt*alpha_adapt*(1*(-a_i[order_i[i]][t]+(u_i[order_i[i]][t]-I_i)));
                        
                        //Spiking activity of inhibitory neurons
                        
                        //Poisson neurons
                        
                        long seed_q =rand()%101;
                        long seed_i = (long) 3+seed_q*i+11*i+i+t*q+s1+67*s2;
                        double p_i =ran2(&seed_i);
                        double p_fire_i = ( 1-exp(-f(u_i[order_i[i]][t+1],h_i[order_i[i]])*dt));
                        if (p_i<p_fire_i)//spike occurs
                            
                            
                        {
                            X_i[order_i[i]][t+1] =1/dt;
                            
                            
                            
                            
                        }
                        
                        else//no spike
                        {
                            X_i[order_i[i]][t+1]=0;
                            
                        }
                        
                        
                        
                    }
                    
                    
                    for (int i=0;i<Ne;i++)
                    {
                        mean_u_e[t+1] = mean_u_e[t+1]+1/((double)Ne)*u_e[order_e[i]][t+1];
                        
                    }
                    for (int i=0;i<Ni;i++)
                    {
                        mean_u_i[t+1] = mean_u_i[t+1]+1/((double)Ni)*u_i[order_i[i]][t+1];
                    }
                    
                    /*
                     //EQUIVALENT MEAN FIELD
                     U_e[t+1]=U_e[t]+dt*alpha_e*(-U_e[t]+b*V_e[t]+w_ee*F(U_e[t], sigma_e)+w_ie*F(U_i[t], sigma_i)+I_e)+sqrt(2*(D)*dt)*Xi_e[0][t]/((double)Ne);
                     U_i[t+1]=U_i[t]+dt*alpha_i*(-U_i[t]+b*V_i[t]+w_ei*F(U_e[t], sigma_e)+w_ii*F(U_i[t], sigma_i)+I_i)+sqrt(2*(D)*dt)*Xi_i[0][t]/((double)Ni);
                     V_e[t+1]=V_e[t]+dt*alpha_adapt*(-V_e[t]+U_e[t]-I_e);
                     V_i[t+1]=V_i[t]+dt*alpha_adapt*(-V_i[t]+U_i[t]-I_i);
                     */
                    
                    
                }// t loop
                
                for (int t=0;t<T;t++)
                {
                    for(int n=0;n<Ne;n++)
                    {
                        Spikes_e[n][t]=X_e[n][t];
                        
                    }
                    for (int n=0;n<Ni;n++)
                    {
                    	Spikes_i[n][t]=X_i[n][t];
                    }
                }
                
                
                ////////////////////COMPUTE EEG TRACES//////////////////////////////////////////////////////////////////////////
                
                //no clue
                
                ////////////////////SPECTRAL ANALYSIS//////////////////////////////////////////////////////////////////////////
                /*
                 //fourier transforms and power spectral denstities
                 
                 for (int k=0;k<T;k++)
                 {
                 PSD[k] = mean_u_e[k];
                 }
                 
                 unsigned long nn=T/2;
                 four1(PSD-1, nn,1);
                 
                 
                 for (int k=0;k<T;k++)
                 {
                 PSD[k] = 1/((double)T*T)*(fabs(PSD[k])*fabs(PSD[k])+fabs(PSD[(int)T-k])*fabs(PSD[(int)T-k]));
                 }
                 
                 
                 
                 for (int k=0;k<T;k++)
                 {
                 Mean_PSD[k] = Mean_PSD[k]+1/((double) Trials)*PSD[k];
                 }
                 
                 */
                
                ////////////////////RATE CALCULATIONS//////////////////////////////////////////////////////////////////////////
                /*
                 for (int ne=0;ne<Ne;ne++)
                 {
                 double rate_e_counter=0;
                 for (int t=0;t<T;t++)
                 {
                 if(X_e[ne][t]>0.1)
                 {
                 rate_e_counter=rate_e_counter+1/((double)T*Dt);///((double) Ne);
                 }
                 else{}
                 }
                 
                 }
                 for (int ni=0;ni<Ni;ni++)
                 {
                 double rate_i_counter=0;
                 for (int t=0;t<T;t++)
                 {
                 
                 if(X_i[ni][t]>0.1)
                 {
                 rate_i_counter=rate_i_counter+1/((double) T*Dt);///((double) Ni);
                 }
                 else{}
                 
                 }
                 }
                 double mean_rate_e=0;
                 double mean_rate_i=0;
                 
                 
                 
                 
                 
                 int twindow=20;
                 int range_e=Ne;
                 int range_i=Ni;
                 
                 for (int t=twindow;t<T;t++)
                 {
                 double rate_e_counter=0;
                 double rate_i_counter=0;
                 double rate_fdbck_counter=0;
                 double countse=0;
                 double countsi=0;
                 double countsfdbck=0;
                 
                 for (int time=0;time<twindow;time++)
                 {
                 for (int index=0;index<range_e;index++)
                 {
                 if(X_e[(int)Ne/2+index-(int)range_e/2][t-(int)twindow/2+time]>0.1)
                 {
                 rate_e_counter=rate_e_counter+1;///((double) Ne);
                 countse = countse+1;
                 }
                 else{}
                 }
                 }
                 
                 
                 for (int time=0;time<twindow;time++)
                 {
                 for (int index=0;index<range_i;index++)
                 {
                 if(X_i[(int) Ni/2+index-(int)range_i/2][t-(int)twindow/2+time]>0.1)
                 {
                 rate_i_counter=rate_i_counter+1;///((double) Ni);
                 countsi = countsi+1;
                 }
                 else{}
                 }
                 }
                 
                 rate_e[t]=rate_e_counter/((double) twindow*Dt)*1/((double)range_e);
                 rate_i[t]=rate_i_counter/((double) twindow*Dt)*1/((double)range_i);
                 
                 
                 }
                 
                 */
                
                ////////////Spike coherence/////////////////////////////////////////////
                int Pairings = Ne;
                double X1[T];
                double X2[T];
                double mean_corr =0;
                for (int p=0;p<Pairings;p++)
                {
                    long seed1 = 2342*p+56*q+3;
                    long seed2 = 45*p+556*q+45;
                    int x1 = (ran2(&seed1)*Ne);
                    int x2 = (ran2(&seed2)*Ne);
                    for (int t=0;t<T;t++)
                    {
                        X1[t]=X_e[x1][t];
                        X2[t]=X_e[x2][t];
                        
                    }
                    
                    mean_corr = mean_corr+1/((double)Pairings)*coherence(X1,X2,temporal_window);
                }
                COHERENCE[s1][s2] = COHERENCE[s1][s2]+1/((double) Trials)*mean_corr ;
                
                
                for(int n=0;n<Ne;n++)
                {
                    for(int t=0;t<T;t++)
                    {
                        if(Spikes_e[n][t]>0.1)
                        {
                            
                            outfile<<Sigma_e[s1]<<","<<Sigma_i[s2]<<","<<n<<","<<t<<endl;
                        }
                    }
                }
                
                for(int n=0;n<Ni;n++)
                {
                	for(int t=0;t<T;t++)
                	{
                		if(Spikes_i[n][t]>0.1)
                		{

                			outfile3<<Sigma_e[s1]<<","<<Sigma_i[s2]<<","<<n<<","<<t<<endl;
                		}
                	}
                }

                for(int t=0;t<T;t++)
                {
                    outfile2<<Sigma_e[s1]<<","<<Sigma_i[s2]<<","<<t<<","<<mean_u_e[t]<<","<<mean_u_i[t]<<endl;
                }

                
                
                
                
                
            }//-q loop
            
            
            
            /*
             
             
             
             /////////FIND PEAK FREQUENCIES//////////////////////////////////////////////////////////////////////////
             double max_f=0;
             double max_power=0;
             for(int k =2*2*(T)*Dt; k<2*30*(T)*Dt;k++)
             {
             if (	Mean_PSD[k] >max_power)
             {
             max_power = Mean_PSD[k];
             max_f=k*1/((double) 2* (T)*Dt);
             }
             }
             
             
             
             
             */
            
        }//s1-loop
        
    }//-s2-loop
    
    
    ////the following ranks cells by firing rates for raster display purposes//////////////////////////
    /*
     for (int n=0;n<Ne;n++)
     {
     double cumul_spike_count=0;
     for (int t=0;t<T;t++)
     {
     if(X_e[n][t]>0.1)
     {
     cumul_spike_count=cumul_spike_count+1;///((double) Ne);
     }
     else{}
     }
     Indiv_Rates_e[n]=cumul_spike_count/((double)T*Dt);
     Display_indices_e[n]=n;
     }
     */
    /*
     for (int n=0;n<Ne;n++)
     {
     
     for (int i=0;i<Ne-1;i++)
     {
     double temp;
     int temp2;
     if (Indiv_Rates_e[i]>Indiv_Rates_e[i+1])
     {
     temp = Indiv_Rates_e[i+1];
     Indiv_Rates_e[i+1]= Indiv_Rates_e[i];
     Indiv_Rates_e[i] =temp;
     temp2 = Display_indices_e[i+1];
     Display_indices_e[i+1]= Display_indices_e[i];
     Display_indices_e[i] =temp2;
     
     }
     
     }
     
     }
     
     */
    
    
    /////////////OUTPUT//////////////////////////////////////////////////////////////////////////
    
    
    //        ofstream outfile;
    
    for(int k=0;k<(int) T;k++)
    {
        Freq[k] = k*1/((double) 2* (T)*Dt);
    }
    
    
    /*
     
     outfile.open("Heterogeneity -  Firing rates in time.txt", ios::out);
     
     for(int t=0;t<T;t++)
     {
     
     outfile<<t<<"	"<<rate_e[t]<<"	"<<rate_i[t]<<endl;
     
     
     }
     
     outfile.close();
     */
    /*
     outfile.open("Heterogeneity -  Thresholds.txt", ios::out);
     
     for(int n=0;n<Ne;n++)
     {
     
     outfile<<n<<"	"<<h_e[n]<<"	"<<h_i[n]<<endl;
     
     
     }
     
     outfile.close();
     */
    
    
    //        n=sprintf(output3,"Spiking_p%2.0f_S%d_Trials%d.csv",p_connectivity*100,S,Trials);
    //         outfile.open(output3, ios::out);
    //         for(int s1=0;s1<S;s1++)
    //         {
    //             for(int s2=0;s2<S;s2++)
    //             {
    //                 for(int n=0;n<Ne;n++)
    //                 {
    //                     for(int t=0;t<T;t++)
    //                     {
    //                         if(Spikes_e[s1][s2][n][t]>0.1)
    //                         {
    //
    //                                 outfile<<Sigma_e[s1]<<","<<Sigma_i[s2]<<","<<n<<","<<t<<endl;
    //                        }
    //                     }
    //                 }
    //
    //             }
    //
    //
    //
    //         }
    //        outfile.close();
    
    /*
     outfile.open("Heterogeneity -  Spiking data - i cells.txt", ios::out);
     for(int s1=0;s1<S;s1++)
     {
     for(int s2=0;s2<S;s2++)
     {
     for(int n=0;n<Ni;n++)
     {
     for(int t=0;t<T;t++)
     {
     if(Spikes_i[s1][s2][n][t]>0.1)
     {
     
     outfile<<Sigma_e[s1]<<"	"<<Sigma_i[s2]<<"	"<<n<<"	"<<t<<endl;
     }
     }
     }
     
     }
     
     
     
     }
     outfile.close();
     */
    
//    n=sprintf(output1,"MeanSpikeCoherence_p%2.0f_sige%2.0f_sigi%2.0f_Trials%d_Ie%d_gamma%3.0f.csv", p_connectivity*100, sigmae, sigmai, Trials, (int)(I_e*100), gamma*1000);
//    outfile4.open(output1, ios::out);
//    for(int s1=0;s1<S;s1++)
//    {
//    	for(int s2=0;s2<S;s2++)
//    	{
//    		outfile<<Sigma_e[s1]<<","<<Sigma_i[s2]<<","<<COHERENCE[s1][s2]<<endl;
//    	}
//    }
//    outfile4.close();

    /*
     outfile.open("Heterogeneity - Mean PSD.txt", ios::out);
     for(int k=0;k<(int)2*30*(T)*Dt;k++)
     {
     
     outfile<<Freq[k]<<"	"<<Mean_PSD[k]<<endl;
     
     
     }
     outfile.close();
     */
    
    
    
    /*
     outfile.open("Heterogeneity -  Single unit spiking .txt", ios::out);
     for(int t=0;t<T;t++)
     {
     
     outfile<<t<<"	     "<<u_e[12][t]<<"	"<<0.1*X_e[12][t]<<"	"<<u_e[100][t]<<"	"<<0.1*X_e[100][t]<<endl;
     
     }
     outfile.close();
     */
    /*
     outfile.open("Heterogeneity -  Single unit activities.txt", ios::out);
     for(int t=0;t<T;t++)
     {
     outfile<<t<<"	"<<u_e[(int)3*Ne/4][t]<<"	"<<u_e[(int)Ne/4][t]<<"	"<<u_i[(int)3*Ni/4][t]<<"	"<<u_i[(int)Ni/4][t]<<endl;
     
     }
     outfile.close();
     */
    /*
     outfile.open("Heterogeneity -  Single unit Non Linearity.txt", ios::out);
     for(int t=0;t<T;t++)
     {
     outfile<<t<<"	"<<u_e[(int)3*Ne/4][t]<<"	"<<f(u_e[(int)3*Ne/4][t],threshold)<<"	"<<f_prime(u_e[(int)3*Ne/4][t],threshold)<<"	"<<f_2prime(u_e[(int)3*Ne/4][t],threshold)<<endl;
     
     }
     outfile.close();
     */
    
    /*
     outfile.open("Heterogeneity -  E Spikes.txt", ios::out);
     for(int t=0;t<T;t++)
     {
     for (int n=0;n<Ne;n++)
     {
     if(X_e[Display_indices_e[n]][t]>0.1)
     {
     
     outfile<<t<<"	"<<n<<endl;
     }
     }
     
     }
     outfile.close();
     */
    /*
     outfile.close();
     outfile.open("Heterogeneity -  I Spikes.txt", ios::out);
     for(int t=0;t<T;t++)
     {
     for (int n=0;n<Ni;n++)
     {
     if(X_i[n][t]>0)
     {
     
     outfile<<t<<"	"<<n<<endl;
     }
     }
     
     }
     outfile.close();
     */
    
    cout<<"Simulations complete..."<<endl;
    outfile.close();
    outfile2.close();
    outfile3.close();
    outfile4.close();
    outfile5.close();
    
    
    return 0;
}

float ran2(long *idum)
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0) {                             /* initialize */
        if (-(*idum) < 1)                           /* prevent idum == 0 */
            *idum = 1;
        else
            *idum = -(*idum);                         /* make idum positive */
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k*IQ1) - k*IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k*IQ1) - k*IR1;
    if (*idum < 0)
        *idum += IM1;
    k = idum2/IQ2;
    idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
    if (idum2 < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;                                /* avoid endpoint */
    else
        return temp;
}

double F(double u, double sigma)
{	
    double output;
    double v;
    double v_min=-1;
    double v_max=1;
    double V=1000;
    double dv=fabs(v_max-v_min)/((double) V);
    double sum=0;
    for (int i=0;i<V;i++)
    {
        v=v_min+i*dv;
        sum = sum+dv*1/(1+exp(-beta*gamma*(u+v-threshold)))*1/sqrt(2*Pi*sigma*sigma)*exp(-v*v/(2*sigma*sigma));
    }
    output=sum;
    return output;
}

double f(double u,double h)
{
    double output;
    output = 1/(1+exp(-beta*gamma*(u-h)));
    return output;
    
}


void sort(double Array[],int size)
{
    
    for (int q=0;q<size;q++)
    {
        
        for (int i=0;i<size;i++)
        {
            double temp;
            if (Array[i]>Array[i+1])
            {
                temp = Array[i+1];
                Array[i+1]= Array[i];
                Array[i] =temp;
                
            }
            
        }
    }
    
    
}



void shuffle_indices(int Array[], int size)
{
    int temporary;
    int randomNum;
    int last;
    
    for (last = size; last > 1; last=last-1)
    {
        randomNum = (int) floor(rand() % last);
        temporary = Array[randomNum];
        Array[randomNum] = Array[last - 1];
        Array[last - 1] = temporary;
    }
}

void shuffle(double Array[], int size)
{
    double temporary;
    int randomNum;
    int last;
    
    for (last = size; last > 1; last=last-1)
    {
        randomNum = rand() % last;
        temporary = Array[randomNum];
        Array[randomNum] = Array[last - 1];
        Array[last - 1] = temporary;
    }
}


/******************************************************************************/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
 Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
 1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
 if isign is input as -1.  data is a complex array of length nn or, equivalently,
 a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
 checked for!).
 *******************************************************************************/
{
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;
    
    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
        if (j > i) {
            SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
            SWAP(data[j+1],data[i+1]);
        }
        m=nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    
    mmax=2;
    while (n > mmax) { /* Outer loop executed log2 nn times. */
        istep=mmax << 1;
        theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
            for (i=m;i<=n;i+=istep) {
                j=i+mmax; /* This is the Danielson-Lanczos formula. */
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

double coherence(double X1[], double X2[], int window)
{
    int Bins = T/((int) window);
    double Y1[Bins];
    double Y2[Bins];
    //create binned spike train representations
    for (int k=0;k<Bins;k++)
    {
        Y1[k]=0;
        Y2[k]=0;
        for (int s=0;s<window;s++)
        {
            if (X1[k*window+s]>0)
            {
                Y1[k]=1;
            }
            if (X2[k*window+s]>0)
            {
                Y2[k]=1;
            }
            
        }
        
    }
    //compute correlation
    double sumY1Y2=0;
    double sumY1=0;
    double sumY2=0;
    for (int k=0;k<Bins;k++)
    {
        sumY1Y2 = sumY1Y2+Y1[k]*Y2[k];
        sumY1 = sumY1+Y1[k];
        sumY2 = sumY2+Y2[k];
    }
    
    double output;
    if (sumY1>0&&sumY2>0)
    {
        output = sumY1Y2/sqrt(sumY1*sumY2);
    }
    else{ output =0;}
    return output;
    
    
    
}




