#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"


//#ifdef ATLAS  //  1903.06248  and https://www.hepdata.net/record/ins1725190

int main(void)
{
   sortOddParticles(NULL);

   txtList L;
   double MZp=pMass("Z");
   double wZp=pWidth("Z",&L);   
   double BrEe=findBr(L,"e,E"); 
   printf("Z width=%E, branching{e,E}=%E\n",wZp,BrEe);
   numout* cc=getMEcode(1,0,"P,aP->e,E{d,u{D,U","A","","Zp_dudEe");
//printf("cc=%p\n", cc);   
   addCut(1,"M12", MZp-2*wZp, NAN);
   addCut(1,"T(e)",30,  NAN);
   addCut(1,"Y(e)",-2.5,2.5);
   addCut(1,"Y(E)",-2.5,2.5);
//   addHistogram("Y(E)",-2.5,2.5);
//   addHistogram("M(e,E)",500,1500);
//   addHistogram("T(e)",30,600);

   vegasGrid * vegPtr;
   double cs=0;
   int nsub,nin,nout;
   procInfo1(cc,&nsub,&nin,&nout);  
   for( int k=1;k<=nsub;k++)
   {   
      double ti,tsi;
      vegPtr= initMCsession(cc,k,6.5E3, 6.5E3, "PDT:CT10(proton)", "PDT:CT10(proton)");
      if(vegPtr) for(int i=0;i<5;i++)  vegas_int(vegPtr, 100000, 1.5, 4 ,  &ti, &tsi);
      printf("%s,%s->%s,%s   cs=%.2E +/- %.1E [pb] \n", cc->interface->pinf(k,1,NULL,NULL),
                                                        cc->interface->pinf(k,2,NULL,NULL),
                                                        cc->interface->pinf(k,3,NULL,NULL),
                                                        cc->interface->pinf(k,4,NULL,NULL), 
                                                      ti,tsi);
      cs+=ti;   
      vegas_finish(vegPtr);    
   }        
   printf("MZ=%E\n", MZp);  
//   cs*=2*ZpKNNLO(MZp)*1000;  //  fb
   cs*=1000; 
   printf("cross section = %.2E [fb], number of expected Zp  events = %.2E  \n", 
          cs,cs*139); 
}     
//#endif
