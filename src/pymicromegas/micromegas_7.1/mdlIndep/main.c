#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"



int main(void) 
{  
   sortOddParticles(NULL);
//   numout*cc=newProcess("b,B->A,A,A");
   numout*       cc4A=getMEcode(0,0,"t,T->Z,Z,Z","h",NULL,"tt3Z");
   int err=0;
   double Pcm=0.01;
 //  double cs= cs23(cc, 1, Pcm, 4 ,&err);
   double T=0.001;
   double cs=vSigmaCC(T,cc4A,0);
   printf("vcs=%E sigma~%E \n",cs, cs/(2*Pcm/pMass("t")));
   cs=cs23(cc4A, 1, Pcm,  3,&err);
   printf("cs(cs23=%E\n",cs);
   double dcs,chi2;
   cs=cs23Vegas(cc4A,1, Pcm, 3, 0,5, 100000,  5, 10000, 
    &dcs, &chi2); 
   printf("cs(vegas)=%E\n",cs); 
   
   
}    