#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"

//#define FIG5   

double Mtest[7] = { 8,10,30,100,150,300,1000};

double lim95(double M,  double * vs) { return polint3(M,7,Mtest,vs); } 

int  main(void)
{
  SpectraFlag=2;   

  double vSigma=2.9979E-26; //  [cm^3/s] = 1pbc
  double  SpA[NZ];
#ifdef FIG5               //  reconstruction of exclusion plot  Fig.5/Right 2002.01229
   double vs[7];
   for( int i=0;i<7;i++)
   { 
     basicSpectra(Mtest[i],5,0,SpA);     
     Mcdm=Mtest[i];
     vs[i]=vSigma/DwarfSignal(vSigma,SpA); // depends implicitely on Mdcm
   }
   displayPlot("95%% excluded vSigma","M[GeV]",Mtest[0],Mtest[6],1,1,"vSigma[cm^3/s]",0,lim95,vs); 
    
#else  // single point 

vSigma/=3.586373E+00;
vSigma*=1E-5;
  Mcdm=8;
  basicSpectra(Mcdm,5,0,SpA);  
  double res=DwarfSignal(vSigma,SpA);
  printf("res=%E\n",res);
  printf("DwarfSignal> 1  indicates that the model is excluded at the  95%% level\n");
#endif


  
}