
#include"MOinterface.h"
#include"subproc.h"
#include"rw_sess.h"
#include"kininpt.h"
#include"alphas2.h"
#include"runVegas.h"
#include"dynamic_cs.h"

#include"strfun.h"
#include"q_kin.h"
#include "nType.h"

#include"file_scr.h"
#include"cut.h"
#include"histogram.h";

void  delCuts(void)      {  cleartab(&cutTab);}

void delHistograms(void) {  cleartab(&histTab);}

int saveHist(char * fname)
{
  FILE * d=fopen(fname,"w");  
  wrt_hist2(d,"");
  fclose(d);
}


static  void addLinetoTab(table*tab, char *txt)
{ linerec * lr=malloc(sizeof(linerec));
  strcpy(lr->line,txt);
  lr->pred=NULL;
  lr->next=tab->strings;
  if( tab->strings) 
  {  
     tab->strings->pred=lr;
  }   
  tab->strings=lr;
}   

int  addCut(int key,char*param, double min, double max)
{
   char txt[200];
   if(!isfinite(min) && !isfinite(max)) return 1;
   if(isfinite(min) && isfinite(max) && max<=min ) return 2;
   if( key<0) sprintf(txt,"! |"); else sprintf(txt,"  |");
   sprintf(txt+strlen(txt),"%s |",param);
   if(isfinite(min)) sprintf(txt+strlen(txt),"%E |",min); else sprintf(txt+strlen(txt),"         |");
   if(isfinite(max)) sprintf(txt+strlen(txt),"%E  ",max); else sprintf(txt+strlen(txt),"         ");
   addLinetoTab(&cutTab, txt);
   return 0;
}

int  addHistogram(char*param, double min, double max)
{
   char txt[200];
   if(!isfinite(min) || !isfinite(max)) return 1;
   if(max<=min ) return 2;
   sprintf(txt,"%s |",param);
   sprintf(txt+strlen(txt),"%E |",min);
   sprintf(txt+strlen(txt),"%E |              |          |         | ",max);
   addLinetoTab(&histTab, txt);
   return 0;
}

static char procTxt[100];


int saveHistograms(char *fname)
{  FILE*F=fopen(fname,"w");
   wrt_hist2(F,procTxt);
   fclose(F);
   return 0;
}



vegasGrid * initMCsession( numout*cc, int nSubProc, double P1,double P2, char*STF1, char*STF2)
{ 
  if(cc==NULL) { printf("initMCsession:Squared mutrix element is not compiled\n"); return NULL;} 
  int err=passParameters(cc);
  if(err) {printf("initMCsession:Can not calculate internal constrained parameters for matrix elements\n");  return NULL;}  
  int ntot, nin,nout;  
  procInfo1(cc,&ntot,&nin,&nout);
  if(nin!=2) { printf("initMCsession:We expect a process with 2 incoming particles\n");  return NULL;}
  if(nSubProc<0 || nSubProc>ntot)
             { printf("initMCsession:wrong number of subprocess\n");  return NULL;}
  link_process(cc->interface);
  Nsub=nSubProc;  
  procTxt[0]=0;
  for(int n=1;n<=nin+nout;n++) 
  { sprintf(procTxt+strlen(procTxt),"%s",cc->interface->pinf(Nsub,n,NULL,NULL)); 
    if(n==nin) sprintf(procTxt+strlen(procTxt),"->"); else if(n<nin+nout) sprintf(procTxt+strlen(procTxt),","); 
  }
 
  if(loadStrFun(STF1,STF2)) return NULL; 
  if(fillCutArray()) return NULL;     // indeed the program is terminated    
  if(correctHistList()) return NULL;  // before returns MULL via wrinting a message in the blind mode
  int dim=imkmom(P1,P2);
  clearHists();
  return vegas_init(dim,vegasSQME,50); 
} 
