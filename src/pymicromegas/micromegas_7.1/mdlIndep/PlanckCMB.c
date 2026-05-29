#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

int main(void)
{
   double vcs= 2.9979E-26; // [cm^3/s] = 1[pb c]
   Mcdm = 10; // If model is not loaded,  Mcdm is needsed to calculate the DM number density.
   int pdg =5; // b-quark
   double eSpect[NZ],gSpect[NZ];
   basicSpectra(Mcdm,pdg,0,gSpect);
   basicSpectra(Mcdm,pdg,1,eSpect);
   double exc=PlanckCMB(vcs,gSpect, eSpect);
   printf("Mcdm=%.2e, result=%.2E\n",Mcdm,exc);
   if(exc>1) printf("Excluded\n"); else printf("Not excluded\n");
   return 0;
}
