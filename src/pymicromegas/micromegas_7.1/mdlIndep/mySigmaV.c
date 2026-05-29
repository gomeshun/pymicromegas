#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"





int main(void)
{
  int fast=0; double Beps=1E-4;  // nothing depends on these parameters 
  int err;

  double vs=1; //[pbc]
  sortOddParticles(NULL);

  double Xf=0,T0=1E4;
  Mcdm=400;
  
  constSigmaVpbc=&vs;
  
  double Omega;
  
  Omega=darkOmega(&Xf,fast,Beps,&err); 
  printf("OmegaFO=%.3E Xf=%.3E\n", Omega,Xf);
     
  Omega=darkOmegaTR(T0,0,fast,Beps,&err);
   
  displayPlot("Y","T",Tend,Tstart,1,2,"Y",0,YF,NULL, "Yeq",0,Yeq,NULL);
  
  displayPlot("vs","T",2,Tstart,0,1,"vSigma",0,vSigmaMem,NULL);
  
  printf("err=%d Tstart=%.3E Omega=%.3E\n", err,Tstart,Omega);   
}

