#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

   char * name[5]={"Pythia_6", "Pythia_8", "PPPC", "CosmiXs", "CosmiXs2"};
   char * out[7]={"A","e","p","ne","nm","ml","D"};  

static double Egamma;

     double elInt( double e,void *tab)
     {  double *etab=(double *) tab;
        double me=0.511E-3;
//        if(e<me) return 0;
        if(e<Egamma) return 0;
        double p=sqrt((e+me)*(e+me)-me*me);
        return SpectdNdE(e,etab)* FSRdNdE(Egamma,p,me, 1,1,0);
     }
     

int main(int argc,char** argv)
{ 
   double M=20;       //  DM mass
   int pdg=5;         //  b-quark  
   double xMin=1E-7;   //  xMin for plots 
   double SpA[NZ], 
          SpE[NZ], 
          FSR[NZ];

SpectraFlag=0;
  printf("M=%.2E, PDG=%d\n",M,pdg);
  basicSpectra(M,pdg,1,SpE);
  basicSpectra(M,pdg,0,SpA);
  FSR[0]=M;
  int k=1;
  for(k=1;k<NZ;k++)  if(SpE[k]==0) FSR[k]=0; else break;
  for(;k<NZ;k++) 
  {   
     Egamma=M*exp(Zi(k));  
     FSR[k]=Egamma*simpson_arg(elInt,SpE,Egamma,M,1E-3,NULL);
  }

  
       
  displayPlot("gamma","E",M*xMin,M,1,2,"orig",0,eSpectdNdE,SpA,"FSR",0,eSpectdNdE,FSR);

  Mcdm=M;
  double vcs= 2.9979E-26; // [cm^3/s] = 1[pb c]
  double exc=PlanckCMB(vcs,SpA, SpE);
  printf("Without FSR: PlanckCMB=%.4E  ",exc); if(exc>1) printf("Excluded\n"); else printf("not Excluded\n");
  double Emax=M;
  for(int i=0;i<NZ;i++) if(FSR[i]<SpA[i]) { FSR[i]=SpA[i]; Emax=M*exp(Zi(i));}
  double Efsr=simpson_arg(eSpectdNdE, FSR, 1E-10,Emax,1E-3,NULL);
  printf("Efsr/(2M)=%E\n",Efsr/M/2);
  
  exc=PlanckCMB(vcs,FSR, SpE);
  printf("With    FSR: PlanckCMB=%.4E  ",exc); if(exc>1) printf("Excluded\n"); else printf("not excluded\n");

  displayPlot("gamma","E",M*xMin,M,1,2,"orig",0,eSpectdNdE,SpA,"FSR",0,eSpectdNdE,FSR);
}                                           

