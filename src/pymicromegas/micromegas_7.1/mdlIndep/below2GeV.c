#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

int main(int argc,char** argv)
{ int err,n,i;
  SpectraFlag=0;

double etaPh[NZ],etaPos[NZ];
  double M=1;
  basicSpectra(M,221,0,etaPh);
  basicSpectra(M,221,1,etaPos);
  
  displayPlot("EdNdE eta->photons","E[GeV]",1E-5,M,1,1,"photon",0, eSpectdNdE,etaPh); 
  displayPlot("EdNdE eta->positron","E[GeV]",1E-5,M,1,1,"phositrons",0, eSpectdNdE,etaPos);  

  double Eph=simpson_arg(eSpectdNdE,etaPh,1E-5,M,1E-3,NULL); 
  double Epos=simpson_arg(eSpectdNdE,etaPos,1E-5,M,1E-3,NULL);
  printf("Energy Concervation (eta) 1 ?= %E\n",  (Eph+2*Epos)/(2*M) );

  

  

  double pi0[NZ],pich[NZ],kl[NZ],ks[NZ],kch[NZ];
  basicSpectra(M,111,0,pi0);
  basicSpectra(M,211,0,pich);
  basicSpectra(M,130,0,kl);
  basicSpectra(M,310,0,ks);
  basicSpectra(M,321,0,kch);
  displayPlot(" EdN/dE (photon)","E",0.001,M ,1,5,"pi0",0, eSpectdNdE,pi0
                                            ,"pi+/-",0, eSpectdNdE,pich
                                            ,"KL",0, eSpectdNdE,kl
                                            ,"Ks",0, eSpectdNdE,ks
                                            ,"K+/-",0, eSpectdNdE,kch
                                          );
 double ksPos[NZ]; 
 basicSpectra(M,310,1,ksPos);
                                          
   Eph=simpson_arg(eSpectdNdE,ks,1E-5,M,1E-3,NULL); 
   Epos=simpson_arg(eSpectdNdE,ksPos,1E-5,M,1E-3,NULL);
   printf("Energy Concervation (Ks) 1 ?= %E\n",  (Eph+2*Epos)/(2*M) );
                                           
                                          
  killPlots();
}
