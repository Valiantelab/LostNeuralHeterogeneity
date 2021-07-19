#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>




#define S 20 //size of parameter space
#define umin0 -1 // min of u "grid"
#define umax0 1 // max of u "grid"
#define K 1000 // grid steps

#define w_ee0  1.6   //e->e synaptic strenght
#define w_ie0  -4.7 //i->e synaptic strenght
#define w_ei0  3.0   //e->i synaptic strenght
#define w_ii0  -.13 //i->i synaptic strenght

#define beta 300 //300//activation function non-linear gain

#define gamma 0.016

#define threshold  0.0//activation function threshold (inflexion point)
//#define I_e0 -0.0//Bias cuurent e cells !ADD CURRENT HERE TO TRIGGER SEIZURE LIKE STATE/GAMMA
#define I_i0 -0.5// Bias current i cells

#define Pi 3.14159

//#define s1 0
//#define s2 0
#define sigmae 4.4
#define sigmai 16.75

#define Imin0 -.25
#define Imax0 .25
#define Isteps 250




double Sigma_e[S];
double Sigma_i[S];
double Ge_1;
double Ge_2;
double Gi_1;
double Gi_2;
double U_e_1;
double U_e_2;
double U_i_1;
double U_i_2;



double F(double u, double sigma)
{
    double v_min;
    double v_max;
    double vmin0;
    double vmax0;
    vmin0=-1;
    vmax0=1;
    v_min=vmin0/gamma;
    v_max=vmax0/gamma;
    
    double output;
    double v;
    //    double v_min=-33.3;
    //    double v_max=33.3;
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








int main ()
{
    double w_ee;
    double w_ei;
    double w_ie;
    double w_ii;
    double I_e;
    double I_i;
    double Imin;
    double Imax;
    double umin;
    double umax;

    
    w_ee=w_ee0/gamma;
    w_ei=w_ei0/gamma;
    w_ie=w_ie0/gamma;
    w_ii=w_ii0/gamma;
    I_i=I_i0/gamma;
    Imin=Imin0/gamma;
    Imax=Imax0/gamma;
    umin=umin0/gamma;
    umax=umax0/gamma;

    
    char f1[125];
    
    sprintf(f1, "WC_StabilityAnalysisNEW_Iemin%d_Iemax%d_K%d_U%d_sigmae%d_sigmai%d.csv", (int)(Imin*100), (int)(Imax*100), K, (int)umax,(int)(sigmae*100),(int)(sigmai*100));
    
    FILE *Output;
    Output = fopen(f1, "wt");
    
    for (int I=0; I<Isteps; I++)
    {
        
        Sigma_e[0]=sigmae;
        Sigma_i[0]=sigmai;
        I_e=Imin+((float)I/(float)Isteps)*(Imax-Imin);
        printf("Ie= %f\n", I_e);

        //            fprintf(Output, "%3f, %3f, ,", Sigma_e[s1], Sigma_i[s2]);
        
        for (int k1=0;k1<(K);k1++)
        {
            
            //                printf("k1=%d\n", k1);
            
            for (int k2=0;k2<(K);k2++)
            {
                
                //                    printf("k1=%d, k2=%d\n", k1, k2);
                
                // Calculate U_e and U_i for this and the nearest grid points
                U_e_1=umin+((float)k1/(float)K)*(umax-umin);
                U_e_2=umin+(((float)k1+(float)1)/(float)K)*(umax-umin);
                U_i_1=umin+((float)k2/(float)K)*(umax-umin);
                U_i_2=umin+(((float)k2+(float)1)/(float)K)*(umax-umin);
                
                //                    printf("%f, %f, %f, %f, %f, \n", U_e_1, U_e_2, U_i_1, U_i_2, ((float)k1/(float)K)*(umax-umin));
                
                
                //EQUIVALENT MEAN FIELD
                //                    Ge=(-U_e+w_ee*F(U_e, sigma_e)+w_ie*F(U_i, sigma_i)+I_e);
                //                    Gi=(-U_i+w_ei*F(U_e, sigma_e)+w_ii*F(U_i, sigma_i)+I_i);
                
                // Move one grid unit in the U_e direction
                Ge_1=(-U_e_1+w_ee*F(U_e_1, Sigma_e[0])+w_ie*F(U_i_1, Sigma_i[0])+I_e);
                Ge_2=(-U_e_2+w_ee*F(U_e_2, Sigma_e[0])+w_ie*F(U_i_1, Sigma_i[0])+I_e);
                Gi_1=(-U_i_1+w_ei*F(U_e_1, Sigma_e[0])+w_ii*F(U_i_1, Sigma_i[0])+I_i);
                Gi_2=(-U_i_1+w_ei*F(U_e_2, Sigma_e[0])+w_ii*F(U_i_1, Sigma_i[0])+I_i);
                //                    printf ("Ge_1=%f, Ge_2=%f \n", Ge_1, Ge_2);
                
                if (Ge_1*Ge_2<0)
                {
                    fprintf(Output, "%3f,%3f,%3f, 0, %3f,%3f,%3f,%3f,\n", Sigma_e[0], Sigma_i[0],I_e, U_e_1, U_i_1, U_e_2, U_i_1);
                    //                        printf("WE FOUND ONE!!!\n");
                    //                        printf ("U_e_1=%f, U_i_1=%f, Ge_1=%f, Ge_2=%f \n", U_e_1, U_i_1, Ge_1, Ge_2);
                    
                }
                if (Gi_1*Gi_2<0)
                {
                    fprintf(Output, "%3f,%3f,%3f, 1, %3f,%3f,%3f,%3f,\n", Sigma_e[0], Sigma_i[0],I_e,U_e_1, U_i_1, U_e_2, U_i_1);
                    //                        printf("WE FOUND ONE!!!\n");
                    //                        printf ("U_e_1=%f, U_i_1=%f, Ge_1=%f, Ge_2=%f \n", U_e_1, U_i_1, Ge_1, Ge_2);
                    
                }
                
                
                // Move one grid unit in the U_i direction
                Ge_1=(-U_e_1+w_ee*F(U_e_1, Sigma_e[0])+w_ie*F(U_i_1, Sigma_i[0])+I_e);
                Ge_2=(-U_e_1+w_ee*F(U_e_1, Sigma_e[0])+w_ie*F(U_i_2, Sigma_i[0])+I_e);
                Gi_1=(-U_i_1+w_ei*F(U_e_1, Sigma_e[0])+w_ii*F(U_i_1, Sigma_i[0])+I_i);
                Gi_2=(-U_i_2+w_ei*F(U_e_1, Sigma_e[0])+w_ii*F(U_i_2, Sigma_i[0])+I_i);
                //                    printf ("Ge_1=%f, Ge_2=%f \n", Ge_1, Ge_2);
                
                if (Ge_1*Ge_2<0)
                {
                    fprintf(Output, "%3f,%3f,%3f, 0, %3f,%3f,%3f,%3f,\n", Sigma_e[0], Sigma_i[0],I_e,U_e_1, U_i_1, U_e_1, U_i_2);
                    //                        printf("WE FOUND ONE!!!\n");
                    //                        printf ("U_e_1=%f, U_i_1=%f, Ge_1=%f, Ge_2=%f \n", U_e_1, U_i_1, Ge_1, Ge_2);
                    
                }
                if (Gi_1*Gi_2<0)
                {
                    fprintf(Output, "%3f,%3f,%3f, 1, %3f,%3f,%3f,%3f,\n", Sigma_e[0], Sigma_i[0],I_e,U_e_1, U_i_1, U_e_2, U_i_1);
                    //                        printf("WE FOUND ONE!!!\n");
                    //                        printf ("U_e_1=%f, U_i_1=%f, Ge_1=%f, Ge_2=%f \n", U_e_1, U_i_1, Ge_1, Ge_2);
                    
                }
                
                
                //                    // Move one grid unit in the U_e and U_i direction
                //                    Ge_1=(-U_e_1+w_ee*F(U_e_1, Sigma_e[s1])+w_ie*F(U_i_1, Sigma_i[s2])+I_e);
                //                    Ge_2=(-U_e_2+w_ee*F(U_e_2, Sigma_e[s1])+w_ie*F(U_i_2, Sigma_i[s2])+I_e);
                //                    printf ("Ge_1=%f, Ge_2=%f \n", Ge_1, Ge_2);
                //
                //                    if (Ge_1*Ge_2<0)
                //                    {
                //                        fprintf(Output, "%3f,%3f,%3f %3f, ,", U_e_1, U_e_2, U_i_1, U_i_2);
                //                    }
                
                
                
            }
        }
        
        fprintf(Output, "\n");
        
    }
    
    
    fclose(Output);
    
    return 0;
}
