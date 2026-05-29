#include <sys/utsname.h>
#include <unistd.h>
#include "micromegas.h"
#include "micromegas_aux.h"


typedef struct { numout *cc; double T,p0,Pin; REAL*mAddr[4];} cs22IntSt1;
static  REAL zero=0;

static double cs_integrand(double cs, void* argVoid)
{ REAL p[16]={0};
  REAL C=cs;
  cs22IntSt1*arg=argVoid;
  
  for(int i=0;i<2;i++) p[4*i]=sqrt(pow(*(arg->mAddr[i]),2)+ pow(arg->Pin,2));
  REAL Pout=decayPcm(p[0]+p[4],*arg->mAddr[2],*arg->mAddr[3]);
  for(int i=2;i<4;i++) p[4*i]=sqrt(pow(*(arg->mAddr[i]),2)+ pow(Pout,2));

  p[3]=arg->Pin;  p[7]=-p[3];
  p[11]=C*Pout; p[15]=-p[11];
  REAL S=Sqrt(1-C*C);
  p[10]=S*Pout; p[14]=-p[10];
  double GG= sqrt(4*M_PI*alphaQCD(p[0]+p[4]));
  int err;
  return arg->cc->interface->sqme(1,GG,p,NULL,&err);
}

static double  u_integrand( double u, void* argVoid)
{
  if( u==1. || u==0. ) return 0.;

  cs22IntSt1 * arg = ( cs22IntSt1 *) argVoid;
  double m0=*arg->mAddr[0], m1=*arg->mAddr[1];
  double smin=m0+m1;
  double smin_=(*arg->mAddr[2])+(*arg->mAddr[3]);
  if(smin_>smin) smin=smin_;

  double z=u*(2-u);
  double T=arg->T;
  double sqrtS=smin - 3*T*Log(z);
  double J=6*T*(1-u)/z;

  REAL Pcm=decayPcm(sqrtS,m0,m1);
  arg->Pin=Pcm;
  double Pout=decayPcm(sqrtS,*arg->mAddr[2],*arg->mAddr[3]);
  double p0=arg->p0;
    
  int err;
  double cspcms=gauss345_arg(cs_integrand,argVoid,-1,1,1E-3,&err)*Pout/(32*M_PI); //cs*pcm*s^2
  
   double E0=sqrt(p0*p0 + m0*m0);
   double D=p0*sqrt(Pcm*Pcm*sqrtS*sqrtS+m0*m0*m1*m1) +E0*sqrtS*Pcm;
   double p1=(m1*m1*p0*p0 - sqrtS*sqrtS*Pcm*Pcm)/D;
   double intE=exp(-sqrt(p1*p1+m1*m1)/T);   
//printf("s/s'=%E \n", sqrtS*sqrtS/(pow(fabs(p1)+p0,2)-pow(p1+p0,2))  );   
   if(m0!=0) 
   { p1=D/(m0*m0); 
     intE-=exp(-sqrt(p1*p1+m1*m1)/T);
     if(intE<0) intE*=-1;
   }  
   return cspcms*J*intE;
}



double C1Prod(numout*cc,double p, double T)
{
   if(cc==NULL) { printf("C1Prod: Wrong process cc\n"); return 1;}
   int nsub,nin,nout;
   procInfo1(cc, &nsub, &nin, &nout);
   if(nsub!=1 || nin!=2 || nout!=2) { printf("C_chiProd: single 1->2 subrocess is expected\n"); return 2;}

   char*pnames[4];
   REAL*mAddr[4];
   procInfo2(cc,1,pnames,NULL);
   int pid[4];     
   for(int i=0;i<4;i++) 
   { pid[i]=pTabPos(pnames[i]); pid[i]=abs(pid[i])-1;
     char *nm=ModelPrtcls[pid[i]].mass;
     if(nm[0]=='0')mAddr[i]=&zero;
     else  mAddr[i]=varAddress(nm);
   }

   cs22IntSt1 arg;
   arg.cc=cc;
   arg.T=T;   
   arg.p0=p;
   for(int i=0;i<4;i++) arg.mAddr[i]=mAddr[i];
          
   double res=simpson_arg( u_integrand, &arg, 0,1,1E-3,NULL);
   double E0=sqrt(p*p+ pow(*mAddr[0],2));
   res*=T*ModelPrtcls[pid[0]].g*ModelPrtcls[pid[1]].g*exp(-E0/T)/(pow(2*M_PI,5)*p*E0);  
// printf("g=%d\n", ModelPrtcls[pid[0]].g);   
   return res;
}
