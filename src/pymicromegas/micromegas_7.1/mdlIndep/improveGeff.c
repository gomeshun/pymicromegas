#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

#define input  "../Data/hgEff/LM.thg"
#define output "LMi.thg"
#define COMPAR

static double *Ttab,*gTab,*hTab;
static int nTtab;
static double hStd(double T) { return polint3(T,nTtab,Ttab,hTab);}
static void  dgdtDiff(double T, double *G, double *dG)
{  
  double dT=0.001*T;
  double g=*G;
  *dG= (hStd(T)+T/3*(hStd(T+dT)-hStd(T-dT))/(2*dT) - g)*4/T; 
} 
static double zero(double T) { double dT=0.001*T;  return gEff(T)-hEff(T) +T/4*(gEff(T+dT) - gEff(T-dT))/(2*dT) -T/3*(hEff(T+dT) - hEff(T-dT))/(2*dT); }


int main(void)
{
   double *tabs[5];
   int nRec,nCol; 
   nTtab=readTable(input,&nCol,tabs);
printf("nTta%d\n",nTtab);   
/*
   if(tabs[0][0]<tabs[0][nTtab-1]) 
   { for(int i=0;i<nTtab/2;i++) for(int j=0;j<3;j++)     
     { double a1=tabs[j][i],a2=tabs[j][nTtab-1-i];
       tabs[j][i]=a2;tabs[j][nTtab-1-i]=a1;
     }   
   }  
*/      
       
   Ttab=tabs[0];
   hTab=tabs[1];
   gTab=tabs[2];
   FILE*fout=fopen(output,"w");
   fprintf(fout,"# effective degrees of freedom\n"
                "# ****** . Improved for entropy conservation\n"
                "# T[GeV]     heff        geff\n");              
   for(int i=nTtab-1;i>=0;i--)
   {             
      double T=Ttab[i];
      double g;
      if(T<=1E-2 || T>=10) {fprintf(fout,"%.7E %.7E %.7E\n", T,tabs[1][i],tabs[2][i]); g=tabs[2][i];} 
      else
      {  int err,n;
         if(T>0.1 && T<0.25) n=4; else n=2;
         for(int k=0;k<n;k++)
         { 
            double T_=Ttab[i+1];
            double T1=T_*pow(T/T_,(double)k/n), T2=T_*pow(T/T_,(k+1.)/n);
            err= odeint(&g,1, T1, T2,1E-6, fabs(T1-T2)/2, dgdtDiff); 
            fprintf(fout,"%.7E %.7E %.7E\n", T2 ,hStd(T2),g);
         }
      }   
   }
   fclose(fout);

#ifdef COMPAR
#define DIM 200
{
   double T1=1E-2,T2=10;
   double zOld[DIM],zNew[DIM],dg[DIM]; 
   
   loadHeffGeff(output);

   for(int i=0;i<DIM;i++)  { double T=T1*pow(T2/T1,(i+0.5)/DIM); zNew[i]=zero(T); dg[i]=gEff(T); }   
     
   loadHeffGeff(input);
   for(int i=0;i<DIM;i++)  { double T=T1*pow(T2/T1,(i+0.5)/DIM); zOld[i]=zero(T); dg[i]-=gEff(T); dg[i]/=gEff(T); }
      
   displayPlot("Entropy conservation:  gEff(T)-hEff(T) +T/4*(gEff'(T) -T/3*hEff'(T)","T[GeV]", T1,T2,1,2,input, DIM, zOld, NULL
                                                         ,output,DIM, zNew, NULL
                                                                               );      
   displayPlot("gEff correction: gEffNew(T)/gEffOld(T)-1","T[GeV]", T1,T2,1,1,"",DIM,dg,NULL); 
}
#undef DIM   
#endif
   
   
}   
