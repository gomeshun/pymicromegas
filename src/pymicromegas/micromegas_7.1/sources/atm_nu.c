// http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/index.html
#include<stdio.h>
#include<stdlib.h>
#include<math.h> 

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

typedef struct{ float Etab[101]; float fluxTab[4][101][20];} nuFluxStr;

static nuFluxStr*nuFlux=NULL;

char * atmNuFile=NULL;
static char* atmNuFile_cp=NULL;

static int readHondaTab(void) //  (char*name, double*Etab, double fluxTab[4][101][20] )
{ 
  if(atmNuFile_cp) 
  { if(atmNuFile==NULL || strcmp(atmNuFile,atmNuFile_cp)==0)  return 0; 
    else 
    { atmNuFile_cp=realloc(atmNuFile_cp, 1+strlen(atmNuFile));
      strcpy(atmNuFile_cp,atmNuFile);
    }
  }    
  else   
  { if(atmNuFile==NULL) { atmNuFile_cp=malloc(50);  strcpy(atmNuFile_cp,"grn-ally-20-01-mtn-solmax.d");}
      else              { atmNuFile_cp=malloc(1+strlen(atmNuFile)); strcpy(atmNuFile_cp,atmNuFile); }  
  }

  FILE* f=fopen(atmNuFile_cp,"r"); 
  if(f) printf("Reading file %s for atm neutrino\n",atmNuFile_cp);
  else 
  { char* fname=malloc(strlen(micrO)+strlen(atmNuFile_cp)+30);
    sprintf(fname,"%s/Data/data_nu/%s",micrO, atmNuFile_cp);
    f=fopen(fname,"r");
    if(f) printf("Reading file %s for atm neutrino\n",fname);
    free(fname);
  }  
  if(!f) { printf("can not find  file %s\n", atmNuFile_cp);
           free(atmNuFile_cp); atmNuFile_cp=NULL; 
           return 1;
         }     
  
  if(!nuFlux) nuFlux=malloc(sizeof(nuFluxStr));
   
  for(int i=0; i<20;i++)
  {  fscanf(f,"%*[^\n]\n");
     fscanf(f,"%*[^\n]\n");
     for(int j=0;j<101;j++) 
     {  int n;
       if(i==0)
       { n=fscanf(f," %f", nuFlux->Etab+j); 
         if(n!=1) { printf("Error readind energy line %d   \n",j+3);  goto fErr;} 
       } else 
       { float E;
         n=fscanf(f," %f",&E);
         if(n!=1) { printf("Error readind energy  line %d  \n",i*103+j+3); goto fErr; }
         if(E!=nuFlux->Etab[j]) {printf("Energy disagreement line  %d (%E %E)  \n",i*103+j+3,E,nuFlux->Etab[j]); goto fErr; }
       }      
       for(int k=0;k<4;k++)   fscanf(f," %f ",&(nuFlux->fluxTab[k][j][i]));  
     }
  }
  fclose(f);
  return 0; 
  fErr: 
    fclose(f);
    free(atmNuFile_cp); atmNuFile_cp=NULL;
    return 2;  
}  


double atmNuFluxes(int nuPdg, double E, double cs)
{  
  if(readHondaTab()) return NAN;

  int k;
  switch(nuPdg)
  { case  14: k=0; break;
    case -14: k=1; break;
    case  12: k=2; break;
    case -12: k=3; break;
    default: return NAN;
  }
  if(E<0.1 || E>1E4) return NAN;

  int j =(100*log(E) - 100*log(0.1)) / ( log(1E4)-log(0.1)); 
  double alpha=(log(E)-log(nuFlux->Etab[j]))/(log(nuFlux->Etab[j+1]) -log(nuFlux->Etab[j]));
  double css[20],f[20];
  for(int i=0;i<20;i++) { f[i]=nuFlux->fluxTab[k][j][i]; css[i]=0.95-i*0.1; }
  double fj= polint3(cs,20,css,f);
  for(int i=0;i<20;i++) { f[i]=nuFlux->fluxTab[k][j+1][i]; }
  double fj1=polint3(cs,20,css,f);  
  return (1-alpha)*fj+alpha*fj1;     
}

double atmNuFluxesI(int nuPdg, double E)
{
  if(readHondaTab()) return NAN;
  int k;
  switch(nuPdg)
  { case  14: k=0; break;
    case -14: k=1; break;
    case  12: k=2; break;
    case -12: k=3; break;
    default: return NAN;
  }
  if(E<0.1 || E>1E4) return NAN;

  int j =(100*log(E) - 100*log(0.1)) / ( log(1E4)-log(0.1)); 
  double alpha=(log(E)-log(nuFlux->Etab[j]))/(log(nuFlux->Etab[j+1]) -log(nuFlux->Etab[j]));

  double sj=0,sj1=0;
  for(int i=0;i<20;i++) { sj+=nuFlux->fluxTab[k][j][i]; sj1+=nuFlux->fluxTab[k][j+1][i];} 
  return 0.5*0.1*((1-alpha)*sj+alpha*sj1);
}

