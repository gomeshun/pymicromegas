#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

///#define tildG
//#define Mh125
//#define Preliminary
#define Finally

///#define hCc
//#define plot

double set2(double s)
{  double X[6]={303,   424,600,786, 967, 1320};
   double Y[6]={637.5, 534,423,342,289.5, 234};
   double x=162+s/6.*(1320-162);
   double y=polint3(x,6,X,Y);
   double gamma=2*(876-y)/(876-237);
   return gamma*7./3.*1E-7; 
}




int nLoop=5;

 
   double gb3[5]= {0.512938, 2.988810, 10.265142, 22.331087, 17.871924}; 
   double gb4[5]= {0. };
   double gb5[5]= {0.372215, 1.769896,  3.940999,  0.485222, -27.595818}; 



double betaQCD(int nf,double alpha)   // 1606.08659 Eq.6 
{
   double be[5];
   double a=alpha/M_PI;
  
   be[0]= 11./4.  - 1./6.*   nf;
   be[1]= 6.375   - 0.791667*nf;
   be[2]= 22.3203 - 4.36892 *nf  + 0.0940394*pow(nf,2);
   be[3]= 114.23  - 27.1339 *nf  + 1.58238  *pow(nf,2) + 0.0058567*pow(nf,3);
   be[4]= 524.558 - 181.799 *nf  + 17.156   *pow(nf,2) - 0.22586  *pow(nf,3) - 0.0017993*pow(nf,4);

   double s=0;
   for(int k=0;k<nLoop;k++) s+=be[k]*pow(a,k+2);
   return -s*M_PI;
}  

double gammaQCD(int nf, double alpha)  //  1402.6611  Eq.4.1
{ 
 double ga[5];
 ga[0]= 1; 
 ga[1]= 4.20833  - 0.138889*nf;
 ga[2]= 19.5156  - 2.28412 *nf - 0.0270062*pow(nf,2);
 ga[3]= 98.9434  - 19.1075 *nf + 0.276163 *pow(nf,2) + 0.00579322*pow(nf,3);
 ga[4]= 559.7069 - 143.6864*nf + 7.4824   *pow(nf,2) + 0.1083*    pow(nf,3) - 0.000085359*pow(nf,4);
 double a=alpha/M_PI;
 double s=0;
 for(int k=0;k<nLoop-1;k++) s+=pow(a,k+1)*ga[k];
 return -s;
}

double zeta_a(int nf, double alpha)                          // hep-ph/0512060   Eq.43
{ double a=alpha/M_PI;
  double nl=nf-1; 
  double C=0;
  switch(nLoop)
  { case 5: C+= (5.17035 - 1.00993*nl - 0.0219784*nl*nl)*pow(a,4);
    case 4: C+= (0.972057 - 0.0846515* nl)*pow(a,3);
    case 3: C+= 0.152778* pow(a,2);
    case 2:
    case 1: C+=1;
  }   
  return C;
}

double  dAlpha(int nf, double alpha)                  // d alpha*1+zeta(alpha)/d alpha  
{ double a=alpha/M_PI;
  double nl=nf-1; 
  double C=0;
  switch(nLoop)
  { case 5: C+= (5.17035 - 1.00993*nl - 0.0219784*nl*nl)*pow(a,4)*5;
    case 4: C+= (0.972057 - 0.0846515* nl)*pow(a,3)*4;
    case 3: C+= 0.152778* pow(a,2)*3;
    case 2:
    case 1: C+=1;
  }   
  return C;
}


double  zeta_m(int nf,double alpha)                              //  hep-ph/09708255  Eq.20  
{
  double a=alpha/M_PI;
  int  nl=nf-1;
  double s=0;
  switch(nLoop)
  { case 5:
    case 4:s+=(1.8476 + 0.0247*nl)*pow(a,3);
    case 3:s+=0.2060*a*a;   
    case 2:
    case 1:s+=1;  
  }
  return s;   
}

double  dzeta_m(int nf,double alpha)                               
{
  double a=alpha/M_PI;
  int  nl=nf-1;
  double s=0;
  switch(nLoop)
  { case 5:
    case 4:s+=3*(1.8476 + 0.0247*nl)*pow(a,2);
    case 3:s+=2*0.2060*a;   
    case 2: 
    case 1:
  }
  return s/M_PI;   
}

double poleMassF(int nf, double alpha)
{ double sumMass[5]={4.7E-3,7E-3,0.1,1.38,5.56};
  int  nl=nf-1;
  double a=alpha/M_PI;
  double s=0;
  switch(nLoop)
  { case 5:
    case 4:
    case 3:s+=(0.6527*nl*nl- 26.655*nl + 190.595)*pow(a,3);
    case 2:s+=(-1.0414*(nl+4./3.*sumMass[nl])+ 13.4434)*pow(a,2);
    case 1:s+=4./3.*a;
    case 0:s+=1;         
  }
  return s; 
} 
static int nf_qcd=5;

static void betaDer(double lgmu, double*Y, double*dY)
{
  double alpha=Y[0];
  double m=Y[1];
  dY[0]=2*betaQCD(nf_qcd,alpha);
  dY[1]=2*gammaQCD(nf_qcd,alpha)*m;
} 

double solveQCD(int nf, double mu1,double alpha1, double mu2, double *alpha2 )
{
   double x1=log(mu1), x2=log(mu2);
   double Y[2];
   Y[0]=alpha1;
   Y[1]=1;
   nf_qcd=nf;
   odeint(Y,2,x1,x2,1E-6,0.01,betaDer);
   *alpha2=Y[0];
   return  Y[1];  
}


double lightH2G(double Mh, int nL, int cm)
{   
   
   nLoop=nL;
   double MZ=91.187, alphaMZ=0.12, mtmt=164,mbmb=4.18,mcmc=1.27;
   
   double mu[7]={0,0,0, 0,mcmc,mbmb,mtmt};
   
   double R[7][2][2];
   double T[7][2][2];
   
   double alpha,alpha_;
   solveQCD(5, MZ, alphaMZ,  mtmt,&alpha_);
   
   alpha=alpha_; for(int k=0;k<4;k++) alpha= alpha_/zeta_a(6,alpha);
   
   for(int n=6;n>=4;n--) 
   {  alpha_=alpha*zeta_a(n,alpha);
      T[n][0][0]=dAlpha(n,alpha);
      T[n][0][1]=2*dzeta_m(n,alpha);  
      T[n][1][0]=( dAlpha(n,alpha)*betaQCD(n,alpha) - betaQCD(n-1,alpha_))/( 1-2*gammaQCD(n,alpha));  
      T[n][1][1]= zeta_m(n, alpha)+ 2*(gammaQCD(n,alpha) -gammaQCD(n-1,alpha_) +betaQCD(n,alpha)*dzeta_m(n,alpha))
                                          / ( 1-2*gammaQCD(n,alpha)); 
//printf("%d alpha=%E  alpha_=%E\n", n,alpha,alpha_);  
      if(n>4)
      {   int nl=n-1;
          solveQCD(nl, mu[n] , alpha_, mu[nl] ,&alpha);
      
          R[nl][0][0]= betaQCD(nl,alpha)/betaQCD(nl,alpha_);
          R[nl][0][1]= 2*(gammaQCD(nl,alpha) - gammaQCD(nl,alpha_))/betaQCD(nl,alpha_);
          R[nl][1][0]=0;                                    
          R[nl][1][1]=1;
      }
      else alpha=alpha_;                                                                                 
   }
   double tC1=0;

   if(cm) T[4][1][0]*=0.75*cabs(HggF( pow(0.5*Mh/1.67,2))); 

   for(int k1=0;k1<2;k1++) for(int k2=0;k2<2;k2++) for(int k3=0;k3<2;k3++) for(int k4=0;k4<2;k4++)
   
   tC1+=T[6][1][k1]*R[5][k1][k2]*T[5][k2][k3]*R[4][k3][k4]*T[4][k4][0];
 
   tC1*=-0.5/betaQCD(3,alpha);   
         
//printf("tC1t=%E tC1b=%E tC1c=%E\n", tC1t,tC1b,tC1c);
//printf("tC1=%E\n",tC1);
   double alpha_h;           
   solveQCD(3, mcmc, alpha, Mh ,&alpha_h);

   double tg3=0;
   for(int i=0; i<nLoop;i++) tg3+=pow(alpha_h,2+i)*gb3[i];
  
   return  pow(tC1,2)*tg3;
}   

double hqqRqcd(int nl, double alpha)   //Baikov, Chetyrkin, Kuhn  arxiv.org/pdf/hep-ph/0511063.pdf
{ double a=alpha/M_PI;
  
  double s=0;
  switch(nLoop)
  { case 4: s+= (39.34 -nl*(220.9 -nl*(9.685 -nl*0.0205)))*pow(a,4);  // printf(" 4 %E\n",(39.34 -nl*(220.9 -nl*(9.685 -nl*0.0205)))*pow(M_PI,-4));
    case 3: s+= (164.14-nl*(25.77-nl*0.259))*pow(a,3);                // printf(" 3 %E\n",(164.14-nl*(25.77-nl*0.259))*pow(M_PI,-3));
    case 2: s+= (35.94-1.36*nl)*pow(a,2);                             // printf(" 2 %E\n",(35.94-1.36*nl)*pow(M_PI,-2));
    case 1: s+=5.67*a;                                                // printf(" 1 %E\n",5.67/M_PI); 
    case 0: s+=1;
  }
//if(nLoop==1) printf("alpha=%E a=%e a*5.67=%E\n",alpha, a, a*5.67);  
  return s;
}

double lightH2Gnl4(double Mh) { return lightH2G(Mh,4,1);}
double lightH2Gnl5(double Mh) { return lightH2G(Mh,5,1);}

double  wlightH2G(double Mh) { return lightH2G(Mh,5,1)*2*sqrt(2)*1.16E-5*pow(Mh,3)/M_PI;}

double  wlightH2Gs5(double s) { double Mh=sqrt(s);  return lightH2G(Mh,5,1)*2*sqrt(2)*1.16E-5*pow(Mh,3)/M_PI;}
double  wlightH2Gs4(double s) { double Mh=sqrt(s);  return lightH2G(Mh,4,1)*2*sqrt(2)*1.16E-5*pow(Mh,3)/M_PI;}
double  wlightH2Gs3(double s) { double Mh=sqrt(s);  return lightH2G(Mh,3,1)*2*sqrt(2)*1.16E-5*pow(Mh,3)/M_PI;}
double  wlightH2Gs2(double s) { double Mh=sqrt(s);  return lightH2G(Mh,2,1)*2*sqrt(2)*1.16E-5*pow(Mh,3)/M_PI;}
double  wlightH2Gs1(double s) { double Mh=sqrt(s);  return lightH2G(Mh,1,1)*2*sqrt(2)*1.16E-5*pow(Mh,3)/M_PI;}


 double ms(double Q) { double alpha2; return 0.093*solveQCD(3, 2, 0.308, Q, &alpha2);}

 double malpha3(double M) { double  alpha;  solveQCD(3, 2, 0.308, M, &alpha); return pow(alpha*M,3); }  
 
double ChiralDat[101]={ 3.274601E-08  , 3.420866E-08  , 3.571454E-08  , 4.499843E-08  , 5.894277E-08  , 6.816995E-08  , 7.606610E-08  , 8.330385E-08  , 9.016015E-08  , 9.678104E-08  , 1.032547E-07  , 1.096395E-07  , 1.159763E-07  , 1.248258E-07  , 1.355684E-07  , 1.445676E-07  , 1.531070E-07  , 1.614405E-07  , 1.696773E-07  , 1.778777E-07  , 1.860805E-07  , 1.943128E-07  , 2.025948E-07  , 2.109422E-07  , 2.193677E-07  , 2.278820E-07  , 2.364941E-07  , 2.452117E-07  , 2.540417E-07  , 2.629902E-07  , 2.720627E-07  , 2.812643E-07  , 2.905995E-07  , 3.000727E-07  , 3.096879E-07  , 3.194487E-07  , 3.293587E-07  , 3.394214E-07  , 3.496398E-07  , 3.600170E-07  , 3.705561E-07  , 3.812597E-07  , 3.921307E-07  , 4.031716E-07  , 4.143852E-07  , 4.257738E-07  , 4.373400E-07  , 4.490861E-07  , 4.610144E-07  , 4.731273E-07  , 4.854271E-07  , 4.979158E-07  , 5.105958E-07  , 5.234692E-07  , 5.365380E-07  , 5.498045E-07  , 5.632706E-07  , 5.769385E-07  , 5.908101E-07  , 6.048875E-07  , 6.191726E-07  , 6.336675E-07  , 6.483741E-07  , 6.632943E-07  , 6.784301E-07  , 6.937835E-07  , 7.093563E-07  , 7.251503E-07  , 7.411677E-07  , 7.574101E-07  , 7.738794E-07  , 7.905777E-07  , 8.075066E-07  , 8.246680E-07  , 8.420639E-07  , 8.596959E-07  , 8.775661E-07  , 8.956761E-07  , 9.140278E-07  , 9.326230E-07  , 9.514636E-07  , 9.705512E-07  , 9.898878E-07  , 1.009475E-06  , 1.029315E-06  , 1.049409E-06  , 1.069760E-06  , 1.090368E-06  , 1.111236E-06  , 1.132365E-06  , 1.153757E-06  , 1.175415E-06  , 1.197339E-06  , 1.219532E-06  , 1.241995E-06  , 1.264730E-06  , 1.287739E-06  , 1.311023E-06  , 1.334585E-06  , 1.358426E-06  ,
1.382547E-06};
 
 
extern double  alphaQCDnf(double Q,int nf);

double alphaQCDnfCh(double Q, char * nfch) { int nf; sscanf(nfch,"%d",&nf); return alphaQCDnf(Q,nf);} 

double alphaQCDRK(double Q, char * nfch ) 
{ int nf; sscanf(nfch,"%d",&nf); double al;  solveQCD(nf, 0.7, alphaQCDnfCh(0.7,nfch), Q, &al); return al; }

extern double alphaQCD_(double);
extern double MQRun(int nf,double Q);
extern double MtEff_(double Q);
//extern double Mq_Run(double mass2,double Q);
extern double McPole, MtPole;
extern double initQCD2(double Q0,double alphaS, double McMc, double MbMb, double MtMt);


double Li2int(double y) { if(y<0.01) return 1;   return -log(1-y)/y;} 

double Li2(double x) { return simpson(Li2int,0,x,1E-3,NULL);} 

double A(double be) 
{  return (1+be*be)*( 4*Li2((1-be)/(1+be)) + 2*Li2(-(1-be)/(1+be))
   -3*log((1+be)/(1-be))*log(2/(1+be)) -2*log((1+be)/(1-be))*log(be)
                    )
                   -3*be*log(4/(1-be*be))-4*be*log(be);
}

double Deltat(double be) 
{ 
return  1 * 4./3.*0.11/M_PI *(  A(be)/be+1/(16*pow(be,3))*(3+34*be*be-13*pow(be,4))*log((1+be)/(1-be))+3/(8*be*be)*(7*be*be-1)
                   );

}
int main(void)
{
/*
  displayPlot("Mt","MH", 200,10000,1,2,"MtEff",0,MtEff,NULL,"MtEff_",0,MtEff_,NULL);
  exit(0);

   displayPlot("Deltat","be",0.01,0.999, 0,1, "",0,Deltat,NULL);
   exit(0); 
*/
//   double t_g3(double alpha) { double s=0; for(int i=0;i<5;i++) s+=pow(alpha,2+i)*gb3[i]; return s;}

//   displayPlot("t_g3","alpha",0.1,0.2, 0,1,"",0,t_g3,NULL);
/*
   double R1s(double xx) 
   {  double x=xx*xx/4;
      double lms=log(x);   return 9./4. +1.5*lms + x*(-6-18*lms +x*(-33+36*lms+ x*(280./9. +46./3.*lms
   +x*(146./3. +45./2.*lms))));
   }

   double R1sR0s(double x) { return R1s(x)/pow(1-x*x,1.5);}
   
   displayPlot("R1sR0s","x",0.1,0.9,0,1,"",0,R1sR0s,NULL);

   exit(0); 

   
   double Mt_Eff(double Q)
   {  if(Q<2*MtPole) return MtPole;
      double Mx=MtPole;
      double x=pow(2*MtPole/Q,2);
      for(int i=0;i<225;i++) 
      { double Mx1=MtPole* sqrt( (pow(1-4*MtPole*MtPole/Q/Q,1.5)+ alphaQCD(Q)/M_PI*R1s(MtPole*MtPole/Q/Q))/ pow(1-4*Mx*Mx/Q/Q,1.5));     
        Mx=0.9*Mx+0.1*Mx1;
        
        printf("Mx=%E R1s=%E  \n",Mx,R1s(MtPole*MtPole/Q/Q));

      }  
      return Mx;
   }
*/   
   
//exit(0); 


/*
double FX(double tau) {return creal(HggF(tau));}

displayPlot("HggF","tau", 0,10,0,1,"",0,FX,NULL);

   double  ffa(double tau) { return cabs(HggF(tau*tau/4)); }
   double  f43(double tau) { return 4./3.;}
   displayPlot("hgg","t=Mh/Mq",0, 25, 0,2,"|A_1/2(t^2/4)|",0, ffa,NULL
                                      ,"4/3",0,f43,NULL);  


exit(0);
*/
/*
for(int nf=3;nf<7;nf++) 
{
  double   be0= 11./4.  - 1./6.*   nf;
  double   be4= 524.558 - 181.799 *nf  + 17.156   *pow(nf,2) - 0.22586  *pow(nf,3) - 0.0017993*pow(nf,4);
  switch(nf)
  { case 3: printf("%E \n",15.6982- 0.11111*be4/be0); break;
    case 4: printf("%E \n",9.1104 - 0.12000 *be4/be0); break;
    case 5: printf("%E \n",2.69277 - 0.13046 *be4/be0); break;
    case 6: printf("%E \n",-3.5130 - 0.14286 *be4/be0); break;
  }  
}  
exit(0);
*/

/*
   initQCD(0.1180, 1.271, 4.183, 172.56);
   
   Mt_Eff(MtPole*4);
   
   displayPlot("Mt_eff","MH", 170,10000,0,2,"Mt_eff",0,Mt_Eff,NULL,"MtEff",0,MtEff,NULL); 
   printf("MtPole=%E\n",MtPole);

exit(0);
*/
/*
printf("Pole masses : %E %E %E\n",McPole,MbPole,MtPole );
printf(" AlphaS(1.28)=%E\n", alphaQCD(1.28));
*/   
//double Mq_(double Q) { return Mq_Run(1,Q);} 
//   double Mq(double Q) { return MqRun(1,Q);}

//displayPlot("Mq","Q",3,5,0,2,"Mq",0,Mq,NULL,"Mq_",0,Mq_,NULL);
//exit(0);
   
/*   
   printf("alpha = %e\n",alphaQCDnf(91.187,5));
   displayPlot("alphaQCD","Q[GeV]",1,1.5, 0,2,"alpha_S",0, alphaQCD, NULL, "alpha_S",0, alphaQCD_, NULL);
*/
/*
  double McRun2(double Q) { return MQRun(4,Q);}
  displayPlot("McRun","Q",1,200 ,0,1,"McRun",0,McRun2,NULL);   

  displayPlot("Mb","Q",3,200,0,2,"MbEff",0,MbEff,NULL,"MbEff_",0, MbEff_,NULL);   
printf("MbPole=%E\n",MbPole);   
*/
//  displayPlot("McRun","Q",1,200 ,0,1,"McRun",0,McRun,NULL); 
//  printf("alphaQCD(200)=%E\n", alphaQCD(200));
//  initQCD2(200, alphaQCD(200)  , 1.271, 4.183, 162.95);
//   displayPlot("McRun(2)","Q",1,200 ,0,1,"McRun",0,McRun,NULL);
  

  double mqp[3],mqm[3], mmq[3]={1.271, 4.183, 162.95};

 initQCD2(200, 0.10578, 1.271, 4.183, 162.95);

 for(int j=0;j<3;j++) mqp[j]=MQRun(4+j,200.);
    printf("1!mqp=");  for(int j=0;j<3;j++) printf(" %E ",  mqp[j]); printf("\n");
  
  double M[3][3];
  
  double Q0=200;
 
  int nf=3;
  double mu=1, Mh=2;
  
  double alL[3];
  initQCD(0.1180, mmq[0],mmq[1],mmq[2]);

  for(int j=0;j<3;j++) mqp[j]=MQRun(4+j,200.);
    printf("2!mqp=");  for(int j=0;j<3;j++) printf(" %E ",  mqp[j]); printf("\n");

//exit(0);
  double GF=1.166E-5; // GeV^{-2}
  double eps1=1.0001;
  eps1=1.00001;
  for(int i=0;i<3;i++)
  { mmq[i]*=eps1;
    initQCD2(200, 0.10578, mmq[0],mmq[1],mmq[2]);
    double alp=alphaQCDnf(mu,nf); 
    for(int j=0;j<3;j++) mqp[j]=MQRun(4+j,Q0);
//    printf("mqp=");  for(int j=0;j<3;j++) printf(" %E ",  mqp[j]); printf("\n");
    mmq[i]/=eps1*eps1;
    initQCD2(200, 0.10578, mmq[0],mmq[1],mmq[2]);
    double alm=alphaQCDnf(mu,nf);
    for(int j=0;j<3;j++) mqm[j]=MQRun(4+j,Q0);
//    printf("mqm=");  for(int j=0;j<3;j++) printf(" %E ",  mqm[j]); printf("\n");
    for(int j=0;j<3;j++) M[j][i]=log(mqp[j]/mqm[j])/log(eps1)/2;
    alL[i]=log(alp/alm)/log(eps1)/2;
//printf("alp=%E alm=%E\n",alp,alm);    
    mmq[i]*=eps1;
  }    

  for(int i=0;i<3;i++)
  { for(int j=0;j<3;j++)  printf(" %E",M[i][j]); 
    printf("\n");
  }
  
printf("alpha ");  for(int i=0;i<3;i++) printf(" %E",alL[i]);printf("\n");
  
  

  double N[3][3];
  N[0][0]=1/M[0][0];
  N[1][1]=1/M[1][1];
  N[2][2]=1/M[2][2];
  N[1][0]=N[2][0]=N[2][1]=0;
  
  N[0][1]=-M[0][1]*N[0][0]*N[1][1];
  
  N[1][2]=-M[1][2]*N[2][2]*N[1][1];
  N[0][2]=-(M[0][1]*N[1][2] + M[0][2]*N[2][2])*N[0][0];

/*  
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) 
  { double s=0; for(int k=0;k<3;k++) s+=M[i][k]*N[k][j];
    printf(" %d %d -> s %E  M= %E N=%E\n",i,j,s, M[i][j], N[i][j]);
  }
*/  
  alL[0]*=3./4.*creal(HggF(pow(Mh/1.67/2,2)));
  alL[1]*=3./4.*creal(HggF(pow(Mh/5./2.,2))); 
  double C1=0;
  for(int i=0;i<3;i++) for(int k=0;k<3;k++) C1+=alL[i]*N[i][k];
  C1*=-0.25;
  printf("C1=%E\n",C1);  
  if(nf==5)  // comparison
  {  double alpha=alphaQCDnf(mu,nf);    
     double be=betaQCD(nf,alpha); 
     double alpha6=alphaQCDnf(mmq[2],6);
     double gam=gammaQCD(5,alpha6);

     double C_1=-0.5*(betaQCD(nf,alpha)/alpha)*
     (  dAlpha(6,alpha6)*betaQCD(6,alpha6)/betaQCD(5,alphaQCDnf(mmq[2],nf))     -1)
       /(1-2*gammaQCD(6,alpha6) ); 
     
       printf("!  C_1=%E C1*alpha/beta=%E   delta(C1)=%E   \n",C_1, C_1*alpha/be,1-C_1/C1   );
  }       


//printf("alphaQCDnf(nf,mu)=%E beta=%E \n",alphaQCDnf(mu,nf), betaQCD(nf,alphaQCDnf(mu,nf)));  
   double tC1=C1*alphaQCDnf(mu,nf)/betaQCD(nf,alphaQCDnf(mu,nf));

  double alpha=alphaQCDnf(Mh,nf);
  double tg=0;
  if(nf==3) for(int i=0;i<5;i++) tg+=pow(alpha,i+2)*gb3[i];
  if(nf==5) for(int i=0;i<5;i++) tg+=pow(alpha,i+2)*gb5[i];  
  printf("nf=%d w=%E\n",nf,pow(tC1,2)*tg*4*GF*pow(Mh,3)/sqrt(2)/M_PI);
  printf("%E\n", 3./4.*creal(HggF(pow(125./175./2.,2))));
//exit(0);  
  
//   exit(0);

//   displayPlot("ms(Q)","Q", 0.685,0.686,0,1,"ms",0,ms,NULL);
//   exit(0);

     displayPlot("w(M)","M", 1,2,0,1,"lappha^3/M^3",0,malpha3,NULL);
//exit(0);
   double MZ=91.187,mtmt=164,mbmb=4.18,mcmc=1.27;
   double alphaMZ=0.12;
   double alpha_h;
#ifdef tildeG   // calculation gbb=beta^2(alpha)*G(alpha)/alpha^4   

  for(int i=0;i<5; i++) { gb5[i]=0; gb3[i]=0;} 
  G5[5]={1, 3.9523478, 6.9555141, -6.851753,  -75.25914},
  G4[5]={1,1,1,1,1},
  G3[5]={1, 4.6950708, 13.472440, 20.66395,  - 15.96239}; 
  for( int nf=3; nf<=5;nf+=2)   
  {  double be[5];
     be[0]= 11./4.  - 1./6.*   nf;
     be[1]= 6.375   - 0.791667*nf;
     be[2]= 22.3203 - 4.36892 *nf  + 0.0940394*pow(nf,2);
     be[3]= 114.23  - 27.1339 *nf  + 1.58238  *pow(nf,2) + 0.0058567*pow(nf,3);
     be[4]= 524.558 - 181.799 *nf  + 17.156   *pow(nf,2) - 0.22586  *pow(nf,3) - 0.0017993*pow(nf,4);
 
     if(nf==5)  for(int i=0;i<5;i++) for(int j=0;j<5;j++) for(int k=0;k<5;k++)  if(i+j+k<5)   gb5[i+j+k]+=G5[i]*be[j]*pow(M_PI,-j-1)*be[k]*pow(M_PI,-k-1);
     if(nf==3)  for(int i=0;i<5;i++) for(int j=0;j<5;j++) for(int k=0;k<5;k++)  if(i+j+k<5)   gb3[i+j+k]+=G3[i]*be[j]*pow(M_PI,-j-1)*be[k]*pow(M_PI,-k-1);
  } 
  printf("gb3: %f %f %f %f %f \n", gb3[0], gb3[1], gb3[2], gb3[3], gb3[4]); 
  printf("gb5: %f %f %f %f %f \n", gb5[0], gb5[1], gb5[2], gb5[3], gb5[4]); 

#endif   

#ifdef Mh125   // calculation h->gg  for Mh=125   
{  Mh=125;  
   nLoop=5;
   double alpha5Mt, alpha6Mt;
   solveQCD(5, MZ, alphaMZ,  mtmt,&alpha5Mt);
   alpha6Mt=alpha5Mt; for(int k=0;k<4;k++) alpha6Mt= alpha5Mt/zeta_a(6,alpha6Mt);
   
   
   solveQCD(5, MZ, alphaMZ, Mh,&alpha_h);
   
   double  resSTD= pow(alpha_h/(12*M_PI),2)*(1 + 5.703052*alpha_h +15.57384*pow(alpha_h,2) +12.5520 *pow(alpha_h,3) - 72.0916*pow(alpha_h,4)); 
   printf("resSTD          =%e\n",resSTD);
 
   double tC1=-0.5*( dAlpha(6,alpha5Mt)*betaQCD(6,alpha6Mt)/betaQCD(5,alpha5Mt) -1)/(1-2*gammaQCD(6,alpha6Mt)); 
   double tg5=0;
   for(int i=0; i<5;i++) tg5+=pow(alpha_h,2+i)*gb5[i];
   printf("resRGE=tC1^2*tg5=%e\n", tC1*tC1*tg5); 
}
#endif
   
//  2 GeV Higgs (preliminary)

#ifdef Preliminary
  printf(" Mh=2GeV. Preliminary\n");
for(nLoop=1; nLoop<=5; nLoop++)
{
//   double Mh=3;
   double mtmt=164;
   double alpha5Mt, alpha6Mt;
   solveQCD(5, MZ, alphaMZ,  mtmt,&alpha5Mt);
 
   alpha6Mt=alpha5Mt; for(int k=0;k<4;k++) alpha6Mt= alpha5Mt/zeta_a(6,alpha6Mt);
         
   double mbmb=4.18;
   double alpha5Mb, alpha4Mb;
   solveQCD(5, MZ, alphaMZ,  mbmb,&alpha5Mb);
   alpha4Mb=alpha5Mb*zeta_a(5,alpha5Mb); 
   
   double mcmc=1.27;
   double alpha4Mc, alpha3Mc;
   solveQCD(4, mbmb, alpha4Mb, mcmc ,&alpha4Mc);
   alpha3Mc=alpha4Mc*zeta_a(4,alpha4Mc); 

//   printf("alpha_t %E %E\n", alpha6Mt,alpha5Mt);
//   printf("alpha_b %E %E\n", alpha5Mb,alpha4Mb);
//   printf("alpha_c %E %E\n", alpha4Mc,alpha3Mc);

     
   double tC1t=-0.5*( dAlpha(6,alpha5Mt)*betaQCD(6,alpha6Mt)/betaQCD(5,alpha5Mt) -1)/(1-2*gammaQCD(6,alpha6Mt))
                    * dAlpha(5,alpha4Mb)*betaQCD(5,alpha5Mb)/betaQCD(4,alpha4Mb)
                    * dAlpha(4,alpha3Mc)*betaQCD(4,alpha4Mc)/betaQCD(3,alpha3Mc);
                   
   double tC1b=-0.5*( dAlpha(5,alpha4Mb)*betaQCD(5,alpha5Mb)/betaQCD(4,alpha4Mb) -1)/(1-2*gammaQCD(5,alpha5Mb))
                    * dAlpha(4,alpha3Mc)*betaQCD(4,alpha4Mc)/betaQCD(3,alpha3Mc);

   double tC1c=-0.5*( dAlpha(4,alpha3Mc)*betaQCD(4,alpha4Mc)/betaQCD(3,alpha3Mc) -1)/(1-2*gammaQCD(4,alpha4Mc))
                    ;              
   solveQCD(3, mcmc, alpha3Mc, Mh ,&alpha_h);

   double tg3=0;
   for(int i=0; i<nLoop;i++) tg3+=pow(alpha_h,2+i)*gb3[i];

//   printf("tC1t=%E tC1b=%E tC1c=%E\n", tC1t,tC1b,tC1c);
//   printf("tg3=%E\n", tg3);   
   printf(" alpha(Mh)=%E  res_0=%.2E    res_%d-loop=%.2E\n", alpha_h, pow(3*alpha_h/(12*M_PI),2),nLoop,pow(tC1t+tC1b+tC1c,2)*tg3);
 }
#endif 




#ifdef Finally 
{
// double Mh=3; 

  printf(" Mh=%.2EGeV\n",Mh);
  printf("h->gg\n"); 
for(nLoop=1; nLoop<=5; nLoop++)
{
   double tw=lightH2G(Mh,nLoop,1);   
   printf("nL=%d ~w=%.2E\n",nLoop,tw*2*sqrt(2)*1.16E-5*pow(Mh,3)/M_PI);
}
  double alpha=0.308,ms=0.0935;
  if(Mh!=2) { nLoop=5;  ms=ms*solveQCD(3,2,alpha,Mh,&alpha); }
  printf("h->s,s   alpha=%.3F  ms=%f \n",alpha,ms);
  
  for(nLoop=0; nLoop<=4; nLoop++)
  {  
    double R=hqqRqcd(3,alpha);
//    printf("R=%E\n", R);   
    printf("nL=%d  w=%.2E\n",nLoop, 1.16E-5*Mh*ms*ms/(4*sqrt(2)*M_PI)*R);
  }

double whSs(double Mh)
{  double alpha=0.308,ms=0.0935;
   nLoop=5;
   if(Mh!=2) { nLoop=5;  ms=ms*solveQCD(3,2,alpha,Mh,&alpha); }
   nLoop=4;
   double R=hqqRqcd(3,alpha);
   return 1.16E-5*Mh*ms*ms/(4*sqrt(2)*M_PI)*R;
}
 
 double ggPss(double Mh) 
 { return wlightH2G(Mh) + whSs(Mh);}
 
 displayPlot("w[GeV]","Mh",1,3,0,2,"h->gg",0,wlightH2G,NULL,"h->Ss",0,whSs,NULL);  
 
 displayPlot("w[GeV]","Mh",1,2.5,0,2,"h->gg+h->Ss",0,ggPss,NULL,"Chiral Model", 101,ChiralDat,NULL);
 
 
 printf("gg %.3E  ss %.3E\n"
   , simpson(wlightH2G,1.5,2.25,1e-3,NULL)/(2.25-1.5)
   ,  simpson(whSs,1.5,2.25,1e-3,NULL)/(2.25-1.5)
   );

} 
#endif 

#ifdef Plot 

   double  ffr(double tau) { return creal(HggF(tau*tau/4));}
   double  ffi(double tau) { return cimag(HggF(tau*tau/4));}
   double  ffa(double tau) { return cabs(HggF(tau*tau/4)); }
   double  f43(double tau) { return 4./3.;}
   displayPlot("hgg","t=Mh/Mq",0, 2, 0,2,"|A_1/2(t^2/4)|",0, ffa,NULL
                                      ,"4/3",0,f43,NULL);  
#endif   

#ifdef hCc

double mh[7]={1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
double wh[7]={1.5106E-29, 1.3665E-25, 5.8639E-24, 6.7039E-23, 4.3009E-22, 2.0876E-21, 9.1898E-21}; 
double Br=0.111/0.666;

//EE*Mc/(2*MW*SW) = Mc*0.31333/(2*0.474*80.385)=Mc*4.112E-3

double lambda=4.112E-3;


double tw[7];

for(int i=0;i<7;i++) 
{   tw[i]=wh[i]*M_PI/(2*lambda*lambda*pow(mh[i],3))*2/Br;
    printf("Mh=%E tw=%E\n", mh[i],tw[i]);
}     

#endif 

#ifdef plot
displayPlot("Gamma_pi&K[GeV]","s[GeV^2]",1,6,0,6,"7/3 set2",0,set2,NULL,
                           "h->gg(L1)",0,wlightH2Gs1,NULL,
                           "h->gg(L2)",0,wlightH2Gs2,NULL,
                           "h->gg(L3)",0,wlightH2Gs3,NULL, 
                           "h->gg(L4)",0,wlightH2Gs4,NULL,  
                           "h->gg(L5)",0,wlightH2Gs5,NULL);

#endif

                                   
return 0;
}
