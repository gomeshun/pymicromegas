#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

int main(void)
{ 
  double inputTab[11][2]=
  { 
    {0.800E-11,13.0E-5},
    {1.104E-11,10.8E-5},
    {1.523E-11,9.60E-5},
    {2.101E-11,8.85E-5},
    {2.899E-11,8.35E-5},
    {4.000E-11,8.00E-5},
    {5.519E-11,7.73E-5},
    {7.615E-11,7.49E-5},
    {10.506E-11,7.24E-5},
    {14.496E-11,7.00E-5},
    {20.000E-11,6.70E-5}
  };

extern int findAsymm(double T, double * L,double * dNnu,double * dNlch);

/*
double L[3]={-1E-3,1E-3,0}, dNnu[3],dNlch[3];

findAsymm(0.1, L,dNnu, dNlch);
printf("dNu=   %E %E %E\n",dNnu[0],dNnu[1],dNnu[2]);
printf("dNlch= %E %E %E\n",dNlch[0],dNlch[1],dNlch[2]);


double ff( double pi, void * xx)
{ double x= *(double *)xx;
  return pi*pi/(1+exp(pi -x));
}  


double nmu(double x)
{  return simpson_arg(ff, &x,  0, 100, 1E-3,NULL);
} 

displayPlot("n(x)","x",0,4,0,1,"",0,nmu,NULL);
*/
//exit(0);

//     double Omega=darkOmegaNu(0.01, 3, 1, 1E-16 ,0.035, -0.035, 0, 10);

     double Omega=darkOmegaNu(0.027, 10, 1, 4E-14 , 0.035 , -0.035, 0, 10);
     double pT=simpson_arg(fNuSterile,"nus 3", 0,10,1E-3,NULL)/simpson_arg(fNuSterile,"nus 2", 0,10,1E-3,NULL);
     printf("Omega h^2=%.2E <p/T>=%.2E\n",  Omega,pT);
     displayPlot("L/s(T)","T",0.01,3,1,1,"",0,LasymmS,NULL);
     displayPlot("e^2f(e)","e=p/T",0,10,0,1,"",0,fNuSterile,"nus 2");
//exit(0);

/*
  for(int n=7;n<8;n++)
  {  double T0=0.01,Tinit=3;  // [GeV]

     double Omega=darkOmegaNu(T0, Tinit, 2, inputTab[n][0], 0 , inputTab[n][1],0,  7.1);
     double pT=simpson_arg(fNuSterile,"nus 3", 0,10,1E-3,NULL)/simpson_arg(fNuSterile,"nus 2", 0,10,1E-3,NULL);
     printf("sin22=%.3E  L/s(init)=%.3E Omega h^2=%.2E <p/T>=%.2E\n",
     inputTab[n][0], inputTab[n][1], Omega,pT);

     displayPlot("L/s(T)","T",0.01,3,1,1,"",0,LasymmS,NULL);
     displayPlot("e^2f(e)","e=p/T",0,10,0,1,"",0,fNuSterile,"nus 2");    
  }
*/     
  killPlots();
}
