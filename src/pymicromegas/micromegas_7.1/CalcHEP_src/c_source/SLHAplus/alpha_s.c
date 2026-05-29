#include "SLHAplus.h"

/*
   this file contains routines for running QCD coupling, 
   running quark masses,  and  effective quarks masses 
   which  generate correct widths of scalar particles.
*/


// default parameters
       double MtPole=172.56;
static double qMass[7] ={ 0,0,0,      0,       1.273, 4.183, 162.9433};   // Mq(Mq) 

// lambdas,MbPole, McPole and qMass[7] are   obtained by  initQCD(0.1181, 1.273, 4.183, 172.56);

static double lambda[7]={ 8.892479E-02, 0, 0, 3.348549E-01, 2.923281E-01, 2.098690E-01, 8.892479E-02 }; // Lambda[0] is a reserve copy of Lambda[6];
       double MbPole=4.937813;
       double McPole=1.679838;

#define MZ 91.1876


static double qMin=0.7;


static double m_fact(int nf, double alpha1, double alpha2);
static int nfMax=6;

static int nLoop=5;

double betaQCD(double alpha,int nf)   // 1606.08659 Eq.6 
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

double gammaQCD(double alpha, int nf)  //  1402.6611  Eq.4.1
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



// https://pos.sissa.it/archive/conferences/260/010/LL2016_010.pdf
static double alpha5L(double Q, double lambda, int nf)
{ 
  double   b0=  (11-(2./3.)*nf)/4,
           b1=  (102-(38./3.)*nf)/16,
           b2=  (2857./2. -5033./18.*nf +325./54.*nf*nf)/64,
           b3=   114.23-27.1339*nf+1.58238*nf*nf+0.0058567*nf*nf*nf,
           b4=  524.56 -181.8*nf+17.16*nf*nf -  0.22586*nf*nf*nf-0.0017993*nf*nf*nf*nf; // 1606.08659
//  b_i  RG coefficients for  a=alpha/pi  b_i=beta_i*pi^{i+1}             
  double t=2*log(Q/lambda), Lt=log(t), b0t=b0*t;
  double d1=b1/(b0*b0t), d2=b2/(b0*b0t*b0t), d3=b3/(b0*b0t*b0t*b0t);

// http://pdg.lbl.gov/2023/reviews/rpp2022-rev-qcd.pdf
  double res=  (M_PI/b0t)*(1-d1*Lt+ d1*d1*( Lt*Lt-Lt-1) +d2 
   - d1*d1*d1*(Lt*Lt*Lt-2.5*Lt*Lt-2*Lt+0.5) -3*d1*d2*Lt +0.5*d3   
   +( 18*b0*b1*b1*b2*(2*Lt*Lt-Lt-1)+ pow(b1,4)*(6*pow(Lt,4)-26*pow(Lt,3)-9*Lt*Lt+24*Lt+7) - b0*b0*b3*b1*(12*Lt+1)+2*b0*b0*(5*b2*b2+b0*b4))/(6*pow(b0*b0t,4))
   );
   
  if(!isfinite(res)) return 0;
  else return res;
}

static double zeta_a(double alpha, int nf)                          // hep-ph/0512060   Eq.43
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

static double  zeta_m(double alpha, int nf)                              //  hep-ph/09708255  Eq.20  
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


static double findLambda(int nf,double alpha, double M)
{ double l1=0.1, l2=0.3;
  double l,a,a1,a2;

  l2=M*exp(-2*M_PI/alpha/(11-2*nf/3.));
  while((a2=alpha5L(M,l2,nf)-alpha) < 0) l2*=1.2;
  l1=l2*0.8;
  while((a1=alpha5L(M,l1,nf)-alpha) > 0)  l1*=0.8;

  do{ l=(l1*a2-l2*a1)/(a2-a1);  
      a=alpha5L(M,l,nf)-alpha;
      if(a<0) { a1=a;l1=l;} else {a2=a;l2=l;}
    } while (fabs(a) > 0.0000001*alpha);
  return l;
}      


static double poleQmass(double mm, double alpha, int nf)
{
  double 
         zeta2=1.6439,
         c1=4./(3.),
         c2=13.4434-1.0414*(nf-1),
         c3=190.595-(nf-1)*(26.655-(nf-1)*0.6527);
  double a=alpha/M_PI;
  switch(nf)
  { case 4: c3=0; break;
    case 5: c2+=zeta2*(1.61/4.62); break; /* (Mc+Ms)/Mb contribution */
    case 6: c2+=zeta2*(6.23/175); /* (Mc+Ms+Mb)/Mt contribution */
  }
  return mm*(1+a*(c1 +a*(c2 +a*c3)));   

}


int  initQCD(double MZalphaS, double McMc, double MbMb, double MtP)
{ 
  lambda[5]=findLambda(5,MZalphaS, MZ);
  MtPole=MtP;

  double MtMt=MtP;
  double alpha6=alpha5L(MtMt,lambda[5],5);

  for(int i=0;i<5;i++)
  { 
    alpha6=alpha5L(MtMt,lambda[5],5)/zeta_a(alpha6,6); 
    MtMt*=MtP/poleQmass(MtMt, alpha6,6) ;
  }
  qMass[6]= MtMt;
//printf("MtMt=%E\n", qMass[6]);  
  lambda[6]= findLambda(6,alpha6,qMass[6]);
  
  nfMax=6;

  double alpha;  
  qMass[5]=0; qMass[4]=0;
  if(MbMb<=lambda[5]) { return lambda[5]; }
  
  qMass[5]=MbMb;
  MbPole=poleQmass(MbMb, alpha5L(MbMb,lambda[5] ,5),5);
  alpha=alpha5L(qMass[5],lambda[5],5);
  alpha*=zeta_a(alpha,5); 
  lambda[4]= findLambda(4,alpha,qMass[5]);

  if(McMc<=lambda[4]) { return lambda[4]; }
  McPole=poleQmass(McMc, alpha5L(McMc,lambda[4],4),4);
   
  qMass[4]=McMc;
  alpha=alpha5L(qMass[4],lambda[4],4);  
  alpha*=zeta_a(alpha,4);
  lambda[3]=findLambda(3,alpha,qMass[4]);
  lambda[3];
//   printf("Lambdas 3- 6 %E %E %E %E\n", lambda[3],  lambda[4], lambda[5],  lambda[6]);  
  lambda[0]=lambda[6];
  return 0;
}



int  initQCD5(double MZalphaS, double McMc, double MbMb, double MtP)
{
   int err= initQCD(MZalphaS,McMc, MbMb,MtP);
   lambda[6]=lambda[5]; 
   nfMax=5;
   return err;   
}

int   initQCDH(double Q0,  double alphaS, double McMc, double MbMb, double MtMt)
{ 
  qMass[6]=MtMt, qMass[5]=MbMb, qMass[4]=McMc;
  double Q=Q0,alpha=alphaS; 
  for(int l=6;l>3;l--)
  { 
    lambda[l]=findLambda(l,alpha,Q);
//printf("!lambda[%d]=%E\n", l, lambda[l]);    
    Q=qMass[l];
    alpha=alpha5L(Q,lambda[l],l);
    alpha*=zeta_a(alpha,l);
  }
  lambda[3]=findLambda(3,alpha,Q);  
  MtPole=poleQmass(MtMt, alpha5L(MtMt,lambda[6] ,6),6);
  MbPole=poleQmass(MbMb, alpha5L(MbMb,lambda[5] ,5),5);
  McPole=poleQmass(MbMb, alpha5L(McMc,lambda[4] ,4),4);
  nfMax=6;
  return 0;
}



static int  NF(double Q)
{ 
         if(Q<qMass[4]) return 3;
   else  if(Q<qMass[5]) return 4; 
   else  if(Q<qMass[6]||nfMax==5) return 5; 
   else                 return 6;
}

double alphaQCD(double Q) 
{ 
  if(Q<qMin) Q=qMin;
  double r= alpha5L(Q,lambda[NF(Q)],NF(Q));
//  printf("Q=%E r=%E  NF=%d lambda[NF(Q)]=%E \n",Q,r,NF(Q),lambda[NF(Q)]);
  return r;
}

double alphaQCDnf(double Q,int nf) { return alpha5L(Q,lambda[nf],nf);}


static double m_fact(int nf, double alpha1, double alpha2)  // 1402.6611 
{  double k,c1,c2,c3,c4;
   switch(nf)
   {
     case 3 : k=4./9;    c1=8.950617E-01;   c2=1.371433E+00;   c3= 1.951696E+00; c4= 9.410804E+00;  break;
     case 4 : k=12./25;  c1=1.014133E+00;   c2=1.389207E+00;   c3= 1.090554E+00; c4= 5.830521E+00;  break;
     case 5 : k=12./23;  c1=1.175488E+00;   c2=1.500707E+00;   c3= 1.724858E-01; c4= 1.664568E+00;  break;
     case 6 : k=4./7;    c1=1.397959E+00;   c2=1.793479E+00;   c3=-6.834291E-01; c4=-3.534417E+00;  break;
   }
   double a1=alpha1/M_PI,a2=alpha2/M_PI;

  { double xx=
           pow(a2,k)*(1+a2*(c1+a2*(c2+a2*(c3+a2*c4))))/
          (pow(a1,k)*(1+a1*(c1+a1*(c2+a1*(c3+a1*c4)))));
    return xx;
  }          
}

double MQRun(double Q,int nf)   // for heavy quarks  
{ 
   if(Q<=qMass[nf]) return qMass[nf];
   if(nf==6 || Q<=qMass[nf+1]) 
   { double alpha1=alpha5L(qMass[nf],lambda[nf],nf), alpha2=alpha5L(Q,lambda[nf],nf); 
     return qMass[nf]*m_fact(nf, alpha1, alpha2);
   } else 
   { 
     double alpha1=alpha5L(qMass[nf],lambda[nf],nf), alpha2=alpha5L(qMass[nf+1],lambda[nf],nf); 
     double Mr=qMass[nf]*m_fact(nf, alpha1, alpha2);
     Mr/=zeta_m(alpha5L(qMass[nf+1],lambda[nf+1],nf+1),nf+1);
     return Mr/qMass[nf+1]*MQRun(Q,nf+1);
   }
}

double MqRun( double mass2GeV, double Q) { return  mass2GeV/MQRun(2.,4)*MQRun( Q,4); }

double McRun(double Q) {  return  MQRun(Q,4); }
double MbRun(double Q) {  return  MQRun(Q,5); }  
double MtRun(double Q) {  return  MQRun(Q,6); } 



static double DeltaQCD(double Q)   // M. Spira, Fortsch. Phys. 46 (1998) 203 [arXiv:hep-ph/9705337]
{
 double  a=alphaQCD(Q)/M_PI;
 int nf=NF(Q);
 double res, res1;

 double r=  a*( 5.67 + a*( (35.94-1.36*nf) + a*(164.14-nf*(25.77-nf*0.259) 
 + a*(39.34 -nf*(220.9 -nf*(9.685 -nf*0.0205)))    //Baikov, Chetyrkin, Kuhn  arxiv.org/pdf/hep-ph/0511063.pdf
 )));
 
 return r;
}


double MqEff(double mass2GeV, double Q) { return MqRun(mass2GeV,Q)*sqrt(1+DeltaQCD(Q)); }

double MQEff(double Q, int nf) 
{ double Mp=0;
  switch(nf)
  { case 6: Mp=MtPole; break;
    case 5: Mp=MbPole; break;
    case 4: Mp=McPole; break;
  }   
  if(Q<2*Mp) return Mp;
  double Meff=MQRun(Q,nf)*sqrt(1+DeltaQCD(Q)); 
  double be2= 1-pow(2*Mp/Q,2);
  return Mp*(1-be2) + Meff*be2;
}  


double MtEff(double Q) { return MQEff(Q,6);} 
double MbEff(double Q) { return MQEff(Q,5);}
double McEff(double Q) { return MQEff(Q,4);} 
 
double nfQCD(double Q) {return NF(Q);}

static double   fiRe(double tau)
{
  double x;
  if(tau<1)
  { x=asin(sqrt(tau));
    return x*x;
  }else if(tau==1) return 0;
  else if(tau>1E10)
  { 
     x=log(4*tau);
     return -0.25*(x*x - M_PI*M_PI);
  }
  else
  {
    x=sqrt(1-1/tau);   
    x=log((1+x)/(1-x));
    return -0.25*(x*x - M_PI*M_PI);
  }
}

static double HggFr(double tau) {  return  2*(tau+(tau-1)*fiRe(tau))/(tau*tau); }


double   hWidthCoeff(double Mh, int nl, double*C1 )
{
// copy QCD initialisation 
   double  McMc=qMass[4],MbMb=qMass[5],MtP=MtPole, alphaMZ=alpha5L(MZ,lambda[5],5);
   int nf5=(lambda[5]==lambda[6]);
   
   if(nf5) lambda[6]=lambda[0];
//   double Q0=MtP+50, alpha0=alphaQCDnf(Q0,6);
       double Q0=200, alpha0=alphaQCDnf(Q0,6);  
   double alL[3],M[3][3];


  double GF=1.166E-5; // GeV^{-2}
  double eps1=1.0001;
  
  for(int i=0;i<3;i++)
  { double alp,alm,mqp[3],mqm[3];
  
    qMass[i+4]*=eps1;
    initQCDH(Q0, alpha0, qMass[4],qMass[5],qMass[6]);
    alp=alphaQCDnf(Mh,nl); 
    for(int j=0;j<3;j++) mqp[j]=MQRun(Q0,4+j);
    qMass[i+4]/=eps1*eps1;
    initQCDH(Q0,alpha0,qMass[4],qMass[5],qMass[6]);
    alm=alphaQCDnf(Mh,nl);
    for(int j=0;j<3;j++) mqm[j]=MQRun(Q0,4+j);
    for(int j=0;j<3;j++) M[j][i]=log(mqp[j]/mqm[j])/log(eps1)/2;
    alL[i]=log(alp/alm)/log(eps1)/2;
    qMass[i+4]*=eps1;
  }    

//  for(int i=0;i<3;i++) printf(" %E",alL[i]);printf("\n");
  
  

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


  for(int k=0;k<3;k++) 
  { C1[k]=0;
    for(int i=0;i<3;i++)  C1[k]+=alL[i]*N[i][k];
  }
  
  initQCD(alphaMZ,McMc,MbMb,MtP); 
  double alpha=alphaQCDnf(Mh,nl);
  
  double gb3[5]= {0.512938, 2.988810, 10.265142, 22.331087,  17.871924}; 
  double gb4[5]= {0.439762, 2.332544,  6.692925,  9.025343, -14.520591}; 
  double gb5[5]= {0.372215, 1.769896,  3.940999,  0.485222, -27.595818}; 

//printf("alpha_h=%E\n", alpha);  
  double tg=0;
  switch(nl)
  { case 3: for(int i=0;i<5;i++) tg+=pow(alpha,i+2)*gb3[i]; break;
    case 4: for(int i=0;i<5;i++) tg+=pow(alpha,i+2)*gb4[i]; break;
    case 5: for(int i=0;i<5;i++) tg+=pow(alpha,i+2)*gb5[i]; break;
  }   
  
  double c=-0.25*alpha/betaQCD(alpha,nl)*sqrt(tg);
  
  C1[0]*=3./4.*HggFr(pow(Mh/McPole/2,2)) *c;
  C1[1]*=3./4.*HggFr(pow(Mh/MbPole/2,2)) *c; 
  C1[2]*=3./4.*HggFr(pow(Mh/MtPole/2,2)) *c; 


  if(nf5) lambda[6]=lambda[5]; 
//printf("C=%E tg=%E GF=%E Mh=%E\n", C1[0]+C1[1]+C1[2],tg,GF,Mh);

  return pow(C1[0]+C1[1]+C1[2],2)*4*GF*pow(Mh,3)/sqrt(2)/M_PI;
}

void   hWidthCoeff2(double Mh, double*C1, double *C2 )
{
   for(int i=0;i<3;i++) {C1[i]=0;C2[i]=0;}
   int nl=5;
   if(Mh<4*qMass[4]) nl=3;
   else if(Mh<4*qMass[5]) nl=4;
   
   hWidthCoeff(Mh, nl, C1 ); 
   if(nl>3) hWidthCoeff(Mh, nl-1, C2 );   
//   printf("Mh=%E nl=%d\n",Mh,nl);
}


double bPoleMass(void) { return MbPole;}
double tPoleMass(void) { return MtPole;}
double cPoleMass(void) { return McPole;}



#ifdef TEST
int main(int n, char **args)
{
  double alphaSMZ,mbp,mbmb,mtp;
  double Q; 
  int nf;
  double McMc=1.4;
  sscanf(args[1],"%lf",&alphaSMZ);
  sscanf(args[2],"%lf",&McMc);
  sscanf(args[3],"%lf",&mbmb);
  sscanf(args[4],"%lf",&mtp);
  sscanf(args[5],"%lf",&Q);

  initQCD(alphaSMZ,McMc,mbmb,mtp);

printf("qMass : %E %E %E %E\n",  qMass[3], qMass[4], qMass[5], qMass[6]); 
printf("lambda: %E %E %E %E\n", lambda[3],lambda[4],lambda[5],lambda[6]);
printf("qMin=%E\n",qMin);

  printf("MbPole=%f\n", MbPole);
  printf("MtMt=%f\n", qMass[6]);

  printf("qmass[6]=%f\n",qMass[6]);

  printf("alphaS(%f)=%f\n",Q,alphaQCD(Q));
  printf("MbRun=%f  MbEff=%f \n",MbRun(Q),MbEff(Q)  );
  printf("MtRun=%f  MtEff=%f \n",MtRun(Q),MtEff(Q)  );
}
#endif
