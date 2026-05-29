
#include <sys/utsname.h>
#include <unistd.h>
#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"
#include "../CalcHEP_src/include/rootDir.h"

#define XSTEP 1.1

static double Yeq_(double T)  // Yeq*exp(Mcdm/T) 
{  double heff;
   double X=Mcdm/T;
   heff=hEff(T);
   return (45/(4*M_PI*M_PI*M_PI*M_PI))*X*X*geffDM(T)*sqrt(M_PI/(2*X))/heff;
}

static double *vsAarr=NULL, *vsSarr=NULL, *vs3arr=NULL, *vs4arr=NULL,    
              *TAarr=NULL,  *TSarr=NULL,  *T3arr=NULL,  *T4arr=NULL;
static int     NaExt=0,     NsExt=0,       N3Ext=0,      N4Ext=0;



static double dYExtNew( double T) 
{ 
  double H=sqrt(8*M_PI/3.*M_PI*M_PI/30.*gEff(T))*T*T/MPlanck;
  double s=2*M_PI*M_PI/45*T*T*T*hEff(T);
  double d=0.001; 
  double dlnYeq_dlnT=(log(Yeq(T*(1+d))) - log(Yeq(T*(1-d))))/(2*d);
  double vSigA=0,vSigS=0,vSig3=0,vSig4=0;
  
  if(vsAarr) vSigA=polint3(T,NaExt,TAarr,vsAarr);
  if(vsSarr) vSigS=polint3(T,NsExt,TSarr,vsSarr);
  if(vs3arr) vSig3=polint3(T,N3Ext,T3arr,vs3arr);
  if(vs4arr) vSig4=exp(polint3(log(T),N4Ext,T4arr,vs4arr));

  double vsEff=(2*vSigA + 0.5*vSigS+0.5*vSig3*exp(-Mcdm/T)+ 2*vSig4*exp(-2*Mcdm/T))/3.8937966E8;  

//printf("H=%E s=%E dlnYeq_dlnT=%E vsEff=%E\n", H,s,dlnYeq_dlnT,vsEff);  
  return H/s*dlnYeq_dlnT/vsEff;
}  




static double *lYtab=NULL;
static double *Ttab=NULL;
static int Ntab=0;



static void XderivLnExt(double s3, double *Y, double *dYdx)
{
  double y=Y[0];
  double yeq, yeq_, sqrt_gStar;
  double T,heff,geff;
  double vSigA=0,vSigS=0,vSig3=0,vSig4=0;
  
  T=T_s3(s3);
  yeq=Yeq(T);
  yeq_=Yeq_(T);

  if(vsAarr) vSigA=polint3(T,NaExt,TAarr,vsAarr);
  if(vsSarr) vSigS=polint3(T,NsExt,TSarr,vsSarr);
  if(vs3arr) vSig3=exp(polint3(log(T),N3Ext,T3arr,vs3arr));
  if(vs4arr) vSig4=exp(polint3(log(T),N4Ext,T4arr,vs4arr));


  heff=hEff(T);
  geff=gEff(T);

  *dYdx=MPlanck
      *pow(2*M_PI*M_PI/45.*heff,2./3.)/sqrt(8*M_PI/3.*M_PI*M_PI/30.*geff)
      *(     vSigA*(y*y-yeq*yeq)
        +0.5*vSigS*(y*y-yeq*y)
        -0.5*vSig3*(y*y*exp(-Mcdm/T)-y*y*y/yeq_)
        -    vSig4*(y*y*exp(-2*Mcdm/T)-y*y*y*y/yeq_/yeq_)
       )/3.8937966E8;
//       printf(" %E\n", yeq/yeq_*exp(Mcdm/T));
        
}




static double (*vsExt)( double T);
static double vsExtLn( double lnT) { return log(vsExt(exp(lnT)));}
 

double darkOmegaExtTR(double Tr,double Yr, double (*vs_a)(double), double (*vs_s)(double) 
                                         , double (*vs_3)(double), double (*vs_4)(double)
                                         , int *err )
{

  double Yt, omega;

  int Nt=25;
  lYtab=realloc(lYtab,sizeof(double)*Nt);
  Ttab=realloc(Ttab,sizeof(double)*Nt);
  Ntab=0;

  Tstart=Tr;
  Yt=Yr;
  double eps=0.001, delta= 0.01;

  if(TAarr) free(TAarr); if(vsAarr) free(vsAarr); TAarr=NULL; vsAarr=NULL;
  if(TSarr) free(TSarr); if(vsSarr) free(vsSarr); TSarr=NULL; vsSarr=NULL;
  if(T3arr) free(T3arr); if(vs3arr) free(vs3arr); T3arr=NULL; vs3arr=NULL;
  if(T4arr) free(T4arr); if(vs4arr) free(vs4arr); T4arr=NULL; vs4arr=NULL;
  
  if(vs_a) buildInterpolation(vs_a,Tend,Tstart, -eps,    delta, &NaExt, &TAarr, &vsAarr);
  if(vs_s) buildInterpolation(vs_s,Tend,Tstart, -eps,    delta, &NsExt, &TSarr, &vsSarr);
  if(vs_3) 
  {  vsExt=vs_3;
     buildInterpolation(vsExtLn,log(Tend),log(Tstart),0.005, delta, &N3Ext, &T3arr, &vs3arr);
  }   
  if(vs_4) 
  {  vsExt=vs_4;
     buildInterpolation(vsExtLn,log(Tend),log(Tstart), 0.005, delta, &N4Ext, &T4arr, &vs4arr);
  }  

     
//  printf("NaExt=%d\n",NaExt);

  Ntab=1;
  Ttab[0]=Tstart;
  lYtab[0]=log(Yt);
  double Tend_=Tstart;

  int checkEq=1;
  for(int i=0; ;i++)
  {
    double s3_t,s3_2,Tbeg;
//    double yeq=Yeq(Mcdm/Xt);

//    if(Tend_<1.E-3 || Yt<fabs(deltaY*1E-5)) break;
    Tbeg=Tend_;
    Tend_/=1.01;
    if(Tend_<Tend) Tend_=Tend;
    s3_t=s3_T(Tbeg);
    s3_2=s3_T(Tend_);
    if(checkEq) 
    { double dY=dYExtNew(Tend_); 
//      printf("dY=%E Yeq=%E\n",dY,Yeq(Tend_));     
      if(fabs(dY)<0.01*Yeq(Tend_))  Yt=Yeq(Tend_)+dY; 
      else  if(odeint(&Yt,1 ,s3_t , s3_2 , 1.E-3, (s3_2-s3_t)/2, &XderivLnExt)){ printf("problem in solving diff. equation\n"); goto theEnd;}
    } else 
    if(odeint(&Yt,1 ,s3_t , s3_2 , 1.E-3, (s3_2-s3_t)/2, &XderivLnExt)){ printf("problem in solving diff. equation\n"); goto theEnd;}  
    if(!isfinite(Yt)||FError)  goto theEnd;
    if(Ntab>=Nt)
    { Nt+=20;
      lYtab=realloc(lYtab,sizeof(double)*Nt);
      Ttab=realloc(Ttab,sizeof(double)*Nt);
    }
    lYtab[Ntab]=log(Yt);
    Ttab[Ntab]=Tend_;
    Ntab++;
    if(Tend_==Tend) break;
    Tbeg=Tend_;
  }
  omega=2.742E8*Mcdm*Yt;

  if(err) *err=0;
  return omega;
  theEnd: 
  if(err) *err=16;
  return NAN;
}

static double dYExt(double s3, double (*vs_a)(double), double (*vs_s)(double))
{ double d, dlnYds3,Yeq0X, sqrt_gStar, vSig,vSig0,vSig1,res;;
  double epsY,alpha,yeq;
  double T,heff,geff;
  T=T_s3(s3);
  yeq=Yeq(T);
  if(yeq<=0) return 10;
  epsY=deltaY/yeq;
  if(vs_a)  vSig0=    vs_a(T)/3.8937966E8; else vSig0=0;
  if(vs_s)  vSig1=0.5*vs_s(T)/3.8937966E8; else vSig1=0;
  vSig=vSig0+vSig1;
  if(vSig <=0) return 10;
  alpha=vSig1/vSig;

  heff=hEff(T);
  geff=gEff(T);
  d=0.001*s3;  dlnYds3=( log(Yeq(T_s3(s3+d)))
                        -log(Yeq(T_s3(s3-d))) )/(2*d);

  res= dlnYds3/(pow(2*M_PI*M_PI/45.*heff,0.66666666)/sqrt(8*M_PI/3.*M_PI*M_PI/30.*geff)
      *vSig*MPlanck*(1-alpha/2)*sqrt(1+epsY*epsY))/Yeq(T);
  res=fabs(res);
  if(res>10) return 10;
  return res;
}



static double darkOmega1Ext(double (*vs_a)(double), double (*vs_s)(double),  double * Xf,double Z1,double dZ1)
{
  double X = *Xf;
  double CCX=(Z1-1)*(Z1+1);
  double dCCX=(Z1-1+dZ1)*(Z1+1+dZ1)-CCX;
  double ddY;
  double dCC1,dCC2,X1,X2;


  ddY=dYExt(s3_T(Mcdm/X),vs_a,vs_s);
  if(FError || ddY==0)  return -1;
  if(fabs(CCX-ddY)<dCCX)
  { *Xf=X;
    return Yeq(Mcdm/X)*sqrt(1+ddY);
  }

  dCC1=dCC2=ddY-CCX; ;X1=X2=X;
  while(dCC2>0)
  {
     X1=X2;
     dCC1=dCC2;
     X2=X2/XSTEP;
     X=X2;
     dCC2=-CCX+dYExt(s3_T(Mcdm/X),vs_a,vs_s);
     if(Mcdm/X>1.E5) return -1;
  }

  while (dCC1<0)
  {
     X2=X1;
     dCC2=dCC1;
     X1=X1*XSTEP;
     X=X1;
     dCC1=-CCX+dYExt(s3_T(Mcdm/X),vs_a,vs_s);
  }
  for(;;)
  { double dCC;
    if(fabs(dCC1)<dCCX)
      {*Xf=X1;  return Yeq(Mcdm/X1)*sqrt(1+CCX+dCC1);}
    if(fabs(dCC2)<dCCX || fabs(X1-X2)<0.0001*X1)
      {*Xf=X2;  return Yeq(Mcdm/X2)*sqrt(1+CCX+dCC2);}
    X=0.5*(X1+X2);
    dCC=-CCX+dYExt(s3_T(Mcdm/X),vs_a,vs_s);
    if(dCC>0) {dCC1=dCC;X1=X;}  else {dCC2=dCC;X2=X;}
  }
}


double darkOmegaExt(double * Xf, double (*vs_a)(double),double (*vs_s)(double) 
                               , double (*vs_3)(double),double (*vs_4)(double) )
{
  double Yt,Xt=25,omega=NAN;
  double Z1=1.1;
  double Zf=2.5;
  int i;
  double Tend_;
  int wimp;
  if(WIMPpos(&wimp)) return NAN;
  double McdmMem=Mcdm;
  Mcdm=McdmN[wimp];
  fracCDM[wimp]=1;
  simpson_err=0;
//  MassCut=4*Mcdm;

  int Nt=25;

  lYtab=realloc(lYtab,sizeof(double)*Nt);
  Ttab=realloc(Ttab,sizeof(double)*Nt);
  Ntab=0;

  Yt=  darkOmega1Ext(vs_a,vs_s,  &Xt, Z1, (Z1-1)/5);

  if(Yt<0||FError) { goto theEnd;}

  if(Yt<fabs(deltaY)*1.E-15)
  {
     if(deltaY>0) dmAsymm=1;  else dmAsymm=-1;
     if(Xf) *Xf=Xt;
     omega=2.742E8*Mcdm*deltaY;
  } else
  {
  Tstart=Mcdm/Xt;

  double eps=0.001, delta= 0.01;

  if(vs_a) buildInterpolation(vs_a,Tend,Tstart, eps, delta, &NaExt, &TAarr, &vsAarr);
  if(vs_s) buildInterpolation(vs_s,Tend,Tstart, eps, delta, &NsExt, &TSarr, &vsSarr);
//  printf("NaExt=%d\n",NaExt);

//  MassCut=MPlanck;

  Ntab=1;
  Ttab[0]=Tstart;
  lYtab[0]=log(Yt);
  Tend_=Tstart;

  for(i=0; ;i++)
  {
    double s3_t,s3_2,Tbeg;
    double yeq=Yeq(Mcdm/Xt);

//    if(Tend_<1.E-3 || Yt<fabs(deltaY*1E-5)) break;
    Tbeg=Tend_;
    Tend_/=1.2;
    if(Tend_<Tend) Tend_=Tend;
    s3_t=s3_T(Tbeg);
    s3_2=s3_T(Tend_);
    if(odeint(&Yt,1 ,s3_t , s3_2 , 1.E-3, (s3_2-s3_t)/2, &XderivLnExt)){ printf("problem in solving diff. equation\n"); goto theEnd;}
    if(!isfinite(Yt)||FError)  goto theEnd;
    if(Ntab>=Nt)
    { Nt+=20;
      lYtab=realloc(lYtab,sizeof(double)*Nt);
      Ttab=realloc(Ttab,sizeof(double)*Nt);
    }
    lYtab[Ntab]=log(Yt);
    Ttab[Ntab]=Tend_;
    Ntab++;
    if(Tend_==Tend) break;
    Tbeg=Tend_;
  }

  if(Xf)
  {  double T1,T2,Y1,Y2,dY2,dY1;
     T1=Ttab[0];
     Y1=exp(lYtab[0]);
     dY1=Zf*Yeq(T1)-Y1;
     *Xf=Mcdm/T1;
     for(i=1;i<Ntab;i++)
     { T2=Ttab[i];
       Y2=exp(lYtab[i]);
       dY2=Zf*Yeq(T2)-Y2;
       if(dY2<0)
       {
         for(;;)
         {  double al,Tx,Yx,dYx,Xx;
            al=dY2/(dY2-dY1);
            Tx=al*T1+(1-al)*T2,  /*Yx=al*Y1+(1-al)*Y2,*/ Yx=exp(polint3(Tx,Ntab,Ttab,lYtab)),   dYx=Zf*Yeq(Tx)-Yx;
            if(fabs(dYx)<0.01*Yx)
            { *Xf=Mcdm/Tx;
              break;
            } else  { if(dYx>0) {T1=Tx,Y1=Yx;}  else {T2=Tx,Y2=Yx;} }
         }
         break;
      }
      else {dY1=dY2; T1=T2; Y1=Y2; *Xf=Mcdm/T2;}
    }
  }


  if(Yt<fabs(deltaY*1E-15))
  {
      if(deltaY>0) dmAsymm=1; else dmAsymm=-1;
      omega= 2.742E8*Mcdm*deltaY;
  }
//  Yi=1/( (Mcdm/Xt)*sqrt(M_PI/45)*MPlanck*aRate(Xt,1,0,NULL,NULL,NULL));
  if(deltaY==0)
  { dmAsymm=0;
    omega= 2.742E8*Mcdm*Yt; /* 2.828-old 2.755-new,2.742 -newnew */
  } else
  {  double a,f,z0,Y0;
     if(Yt<fabs(deltaY*1E-15))
     { if(deltaY>0) dmAsymm=1; else dmAsymm=-1;
       omega= 2.742E8*Mcdm*deltaY;
     }
     a=fabs(deltaY);
     if(Yt<a*1.E-5)  f=Yt*Yt/4/a; else f=(sqrt(Yt*Yt+a*a)-a)/(sqrt(Yt*Yt+a*a)+a);
     z0=sqrt(f)*2*a/(1-f);
     Y0=sqrt(z0*z0+a*a);
     dmAsymm=deltaY/Y0;
     omega=2.742E8*Mcdm*Y0;
  }
  }
  theEnd:
  for(int i=1;i<=Ncdm;i++) fracCDM[i]=0;  fracCDM[wimp]=1;
  Mcdm=McdmMem;
  if(vs_a) { free(TAarr); free(vsAarr); TAarr=NULL; vsAarr=NULL;}
  if(vs_s) { free(TSarr); free(vsSarr); TSarr=NULL; vsSarr=NULL;}
  return omega;
}
