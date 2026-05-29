#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"


void printVar(FILE *f)
{
  int i;
  fprintf(f,"\n# Model parameters:\n");
  for(i=0;i<nModelVars;i++) fprintf(f,"%-6.6s   %f\n", varNames[i], varValues[i]);
}


int assignVal(const char * name, double val)
{
  REAL * a=varAddress(name);
  if(a && a<=varValues+nModelVars )  {*a=val; return 0;} else return 1;
}


int  assignValW(const char*name, double  val)
{ 
  if(assignVal(name,val)==1) { printf(" %s not found\n",  name); return 1;}
  return 0; 
}

/*
int findVal(char * name, double * val)
{
  double * va=varAddress(name); 
  if(va){*val=*va; return 0;} else { printf("problem: findVal(%s)\n",name);   return 1;}  
}

double  findValW(char*name)
{
  double val;

  if(findVal(name,&val))
  {  printf(" name %s not found\n",name);
     return 0;
  } else return val;
}

*/

static int readInfl(FILE*f, char*name)
{
    double val; 
    int n;
    char buff[200];
    fscanf(f,"%[^\n]",buff);
    if(name[0]) n=1+sscanf(buff,"%lf",&val);  else  n=sscanf(buff,"%s %lf",name,&val); 
    if(n==1) 
    {
           if(strcmp(name,"FreezeOut")==0) COSMO=FreezeOut;
      else if(strcmp(name,"FreezeIn")==0)  COSMO=FreezeIn; 
      else if(strcmp(name,"Inflaton")==0) COSMO=Inflaton;
    } else if(n==2)
    {   
             if(strcmp(name,"HInfl")==0)     HInfl=val;
       else  if(strcmp(name,"GammaInfl")==0) GammaInfl=val; 
       else  if(strcmp(name,"BrInfl")==0)    BrInfl=val; 
       else  if(strcmp(name,"HsmInfl")==0)   HsmInfl=val;
       else  if(strcmp(name,"alphaInfl")==0) alphaInfl=val;
       else  if(strcmp(name,"wInfl")==0)     wInfl=val;
    } 
}

int readVar(const char *fname)
{
  double val;
  char name[80];
  int n;
  FILE * f=fopen(fname,"r");
  if(f==NULL) return -1;
  COSMO=FreezeOut;
  for(n=1;;n++)
  { if(fscanf(f,"%s",name)!=1) { n=0; break;}    
    if(name[0]=='#') 
    {
      if(name[1]=='!')  readInfl(f,name+2);  else fscanf(f,"%*[^\n]");
      continue;
    }
    if(fscanf(f,"%lf",&val)!=1) break;
    fscanf(f,"%*[^\n]");
    { 
      int err=assignVal(name,val);
      if(err==1) break;
    }
  }
  fclose(f);
  return n;                                         
}

int readVarSpecial(char *fname, int nVar, char ** names)
{
  int * rdOn;
  double val;
  char name[80];
  int n,i,k;
  FILE*f;

  rdOn=malloc(sizeof(int)*nVar);
  
  for(i=0;i<nVar;i++)rdOn[i]=0;
 
  f=fopen(fname,"r");
  if(f==NULL) return -1;
  COSMO=FreezeOut;
  for(n=1;;n++)
  { if(fscanf(f,"%s",name)!=1) { n=0; break;}
    if(name[0]=='#') 
    {
      if(name[1]=='!')  readInfl(f,name+2);  else fscanf(f,"%*[^\n]");
    }
 
    if(fscanf(f,"%lf",&val)!=1) break;
    fscanf(f,"%*[^\n]");
    { int err;
      for(i=0;i<nVar;i++) if(strcmp(names[i],name)==0) {rdOn[i]=1;break;}
      if(i==nVar) break;
      err=assignVal(name,val);
      if(err==1) break;
      
    }
  }
  fclose(f);
  for(i=0,k=0;i<nVar;i++) if(rdOn[i]==0)
  { if(!k){printf("The following parameters keep default values:\n"); k=1;} 
    { printf("%8.8s=%.4E", names[i],findValW(names[i]));
      if(k==4) {printf("\n");k=1;}else k++;
    }   
  }  
  if(k!=1) printf("\n");
  free(rdOn);
  return n;                                         
}

int  assignval_(char * f_name, double * val, int len)
{  char buff[20];
   fName2c(f_name,buff,len);
   return assignVal(buff, *val);
}
                                                                                
int findval_(char*f_name, double*val, int len)
{
  char c_name[20];
  fName2c(f_name,c_name,len);
  return findVal(c_name,val);
}


void assignvalw_(char* f_name, double * val, int len)
{  char buff[20];
   fName2c(f_name,buff,len);
   assignValW(buff, *val);
}

double findvalw_(char*f_name, int len)  /*FORTRAN*/
{
  char c_name[20];
  fName2c(f_name,c_name,len);
  return findValW(c_name);
}

