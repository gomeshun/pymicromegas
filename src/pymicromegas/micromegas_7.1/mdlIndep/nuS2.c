#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

#define epsMax 8
//#define epsN   200

//#define epsMax 16
#define epsN   500

//#define XexpX 

#define zeta3 1.202  





static double chiIntegrand(double x, void *m_)
{
   if(x==0) return 0;
   double m=*(double *)m_;
   double p=-2*log(x);
   return  0.25*p*p/x/pow(M_PI*cosh(0.5*sqrt(p*p+ m*m)),2) ;          
}

static double nIntegrand( double x, void *m_)
{
   if(x==0) return 0;
   double m=*(double *)m_;
   double p=-2*log(x);
   return  1/(x*M_PI*M_PI)*p*p/(exp(sqrt(p*p+m*m))+1);          
}

static double dndmu(double m) 
{ if(m<0.5) {double c0=8.333255E-02, c03=1.129012E-03; return c0-m*m/0.09*c03;}
  if(m>1.5) return pow(m,1.5)*exp(-m)*sqrt(M_PI/2)/2/M_PI/M_PI*(1+ 1.98/m);
  int err;  double res= simpson_arg(chiIntegrand,&m, 0,1,1E-3,&err); if(err) { printf("line 25 m=%e \n",m); displayPlot("chiI","x",0,1E-25,0,1,"",0, chiIntegrand,&m);    exit(0);} else return res;  
} // dNdmu = 0.5*Eq.(22)

static double np(double m) 
{
 if(m<0.4) {double  c0=9.134449E-02, c03=1.540234E-03; return c0-m*m/0.09*c03;}
 if(m>3) return  pow(m,1.5)*exp(-m)*sqrt(M_PI/2)/2/M_PI/M_PI*(1+ 2.05/m);

 int err;  double res= simpson_arg(nIntegrand,&m, 0,1,1E-3,&err);  if(err)   { printf("line 27 m=%E\n",m); exit(0);} else return res; 
}   // Nmu


static double * ChiTabs[4];
static int powChiTabs=0;

static double Lds[3];

   static void rdChiTabs(void) 
   { int nCol;  
     char*fName=malloc(strlen(micrO)+40);
     sprintf(fName,"%s/Data/ChiTable_alltemp.dat",micrO);
     powChiTabs=readTable(fName,&nCol,ChiTabs);
     free(fName);
   }
     

static    double  ChiBB(double T) { if(!powChiTabs) rdChiTabs(); return polint3(T*1000,powChiTabs,ChiTabs[0],ChiTabs[1]);}
static    double  ChiQQ(double T) { if(!powChiTabs) rdChiTabs(); return polint3(T*1000,powChiTabs,ChiTabs[0],ChiTabs[2]);}
static    double  ChiQB(double T) { if(!powChiTabs) rdChiTabs(); return polint3(T*1000,powChiTabs,ChiTabs[0],ChiTabs[3]);}


/*static*/ int findAsymm(double T, double * L,double * dNnu,double * dNlch)  // T^3 normalized asymmetries   
{   double me=5.11E-4, mm=0.105, ml=1.777;                        // masses of charged  leptons
    double nBng=6E-10,   dB= nBng*2*zeta3/M_PI/M_PI*hEff(T)/hEff(0);     // barion asymmetry normalized on T^3   

    double A[64]={0},C[8]={0};
//  A[i*8+ j]
//  j=0,1,2 - neutrino asymmetry 
//  j=3,4,5 - charger lepton asymmetry
//  j=6 - quark asymmetry  u-d  
//  j=7 - quark asymmetry 2u+d
    double a[6]; 
    a[0]=a[1]=a[2]=2*dndmu(0);  // 2 because of  asymmetry.Particles increas, aParticle decreas.
//printf("T=%E\n",T);    
    a[3]=4*dndmu(me/T); a[4]=4*dndmu(mm/T); a[5]=4*dndmu(ml/T); 
// other solution
    a[0]=a[1]=a[2]=2*np(0);
    a[3]=4*np(me/T); a[4]=4*np(mm/T); a[5]=4*np(ml/T); 


//  Eq.15 :     
    A[0   ]= a[0];   A[3   ]= a[3];                   C[0]=L[0];   // definition of Le
    A[8+1 ]= a[1];   A[8+4 ]= a[4];                   C[1]=L[1];
    A[16+2]= a[2];   A[16+5]= a[5];                   C[2]=L[2]; 
//  Eq. 16    
    A[24  ]= 1;       A[24+3]=-1;                A[24+6]=-1;  // e+nu_e^bar -> d,u^bar 
    A[32+1]= 1;       A[32+4]=-1;                A[32+6]=-1;
    A[40+2]= 1;       A[40+5]=-1;                A[40+6]=-1;

// Eq.18
                                                 A[48+6]=ChiQB(T); A[48+7]=ChiBB(T);  C[7]=dB;
// Eq.19    charge conservation
    A[56+3]=-a[3];   A[56+4]=-a[4];  A[56+5]=-a[5];   A[56+6]=ChiQQ(T); A[56+7]=ChiQB(T); 
/*
double AA[64], CC[8];
for(int i=0;i<8;i++) CC[i]=C[i]; for(int i=0;i<64;i++) AA[i]=A[i]; 
*/
    int err=solveLinEq(8,A, C);
/*
printf(" C[0]=%E   C[2]=%E\n",C[0],C[3]);
printf(" C[1]=%E   C[4]=%E\n",C[1],C[4]);
*/   
    for(int i=0;i<3;i++) dNnu[i]=a[i]*C[i];
    if(dNlch) for(int i=3;i<6;i++) dNlch[i-3]=a[i]*C[i];
/*
double CCC[8];
      
    for(int i=0;i<8;i++) 
    { CCC[i]=0;
      for(int j=0;j<8;j++) CCC[i]+=AA[i*8+j]*C[j];
      printf("CCC[%d]=%E CC=%E\n", i,CCC[i],CC[i]);
    }  
*/    
    return err; 
}

static  double phi(double x) { double c= 2*cosh(sqrt(x)/2); return 1/(c*c);  } 

static double DN1integrand(double x, void* arg)
{ if(x==0) return 0;
  double* arg_=(double*) arg;
  double eps= arg_[0];
  double p=-2*log(x);
//  dp=2 dx/x 
  return 2*p*p/x*phi(eps*eps+p*p) ;
}


double DeltaN_1(double eps)  
{ if(eps>10) return 0;
   double r= simpson_arg(DN1integrand, &eps,0,1,1E-4,NULL)/(M_PI*M_PI);
}



static double DN3integrand(double x, void* arg)
{ if(x==0) return 0;
  double* arg_=(double*) arg;
  double eps= arg_[0];
  double p=-2*log(x);
//  dp=2 dx/x 
  return 2*p*p/x*phi(eps*eps+p*p)*(1-6*phi(eps*eps+p*p)) ;
}

double DeltaN_3(double eps)  { if(eps>10) return 0;return simpson_arg(DN3integrand, &eps,0,1,1E-4,NULL)/(6*M_PI*M_PI);}

  double deltaN_1(double x)
  { double X[12]={0,1,2,3,4,5,6,7,8,9,20,40};
    double Y[12]={-1.7918E+00, -1.9394E+00, -2.3391E+00, -2.9077E+00, -3.5786E+00, -4.3130E+00, -5.0893E+00, -5.8955E+00, -6.7243E+00, -7.5706E+00,  -1.7479E+01,-3.6484E+01};    
    if(x<9) return exp(polint3(x,10,X,Y));
    return  exp(polint1(x,3,X+9,Y+9));
  }

  double deltaN_3(double x)
  { double X[12]={0,1,2,3,4,5,6,7,8,9,20,40};
    double Y[12]={-4.0812E+00, -4.1059E+00, -4.3201E+00, -4.7795E+00, -5.4020E+00, -6.1169E+00, -6.8857E+00, -7.6890E+00, -8.5167E+00, -9.3626E+00,-1.9270E+01, -3.8276E+01};
    if(x<9) return exp(polint3(x,10,X,Y));
    return  exp(polint1(x,3,X+9,Y+9));
  }


double DeltaN_1_0, DeltaN_3_0, DeltaN_1M[3],DeltaN_3M[3];

/*static*/ int findAsymm2(double T, double * L, double * chi)  // T^3 normalized asymmetries   
{   double me=5.11E-4, mm=0.105, ml=1.777;                                   // masses of charged  leptons
    double nBng=6E-10,   dB= nBng*2*zeta3/M_PI/M_PI*hEff(T)/hEff(0);         // barion asymmetry normalized on T^3   

    double A[64]={0},C[8]={0};
//  A[i*8+ j]
//  j=0,1,2 - neutrino asymmetry 
//  j=3,4,5 - charger lepton asymmetry
//  j=6 - quark asymmetry  u-d  
//  j=7 - quark asymmetry 2u+d
    double a[6]; 
    a[0]=a[1]=a[2]=DeltaN_1_0;  
    for(int i=0;i<3;i++) a[3+i]=2*DeltaN_1M[i];

//  Eq.15 :     
    A[0   ]= a[0];   A[3   ]= a[3];                   C[0]=L[0]-DeltaN_3_0*pow(chi[0],3) - 2*DeltaN_3M[0]*pow(chi[3],3) ;   // definition of Le
    A[8+1 ]= a[1];   A[8+4 ]= a[4];                   C[1]=L[1]-DeltaN_3_0*pow(chi[1],3) - 2*DeltaN_3M[1]*pow(chi[4],3);
    A[16+2]= a[2];   A[16+5]= a[5];                   C[2]=L[2]-DeltaN_3_0*pow(chi[2],3) - 2*DeltaN_3M[2]*pow(chi[5],3);
//  Eq. 16    
    A[24  ]= 1;       A[24+3]=-1;                A[24+6]=-1;  // e+nu_e^bar -> d,u^bar 
    A[32+1]= 1;       A[32+4]=-1;                A[32+6]=-1;
    A[40+2]= 1;       A[40+5]=-1;                A[40+6]=-1;

// Eq.18
                                                 A[48+6]=ChiQB(T); A[48+7]=ChiBB(T);       C[6]=dB;
// Eq.19    charge conservation
    A[56+3]=-a[3];   A[56+4]=-a[4];  A[56+5]=-a[5];   A[56+6]=ChiQQ(T); A[56+7]=ChiQB(T); C[7]=DeltaN_3M[0]*pow(chi[3],3)
                                                                                              +DeltaN_3M[1]*pow(chi[4],3)
                                                                                              +DeltaN_3M[2]*pow(chi[5],3);
                                                               //                           C[7]=0;    
    int err=solveLinEq(8,A, C);
//    printf("err=%d\n",err);

    for(int i=0;i<6;i++) chi[i]=C[i];

    return err; 
}




#ifdef DMAS2
/*static*/ int findAsymm2(double T, double * L,double * muN,double * muL)  // T^3 normalized asymmetries   
{   double me=5.11E-4, mm=0.105, ml=1.777;                        // masses of charged  leptons
    double nBng=6E-10,   dB= nBng*2*zeta3/M_PI/M_PI*hEff(T)/hEff(0);     // barion asymmetry normalized on T^3   

    double A[64]={0},C[8]={0};
//  A[i*8+ j]
//  j=0,1,2 - neutrino asymmetry 
//  j=3,4,5 - charger lepton asymmetry
//  j=6 - quark asymmetry  u-d  
//  j=7 - quark asymmetry 2u+d
    double a[6]; 
    a[0]=a[1]=a[2]=2*dndmu(0);  // 2 because of  asymmetry.Particles increas, aParticle decreas.
//printf("T=%E\n",T);    
    a[3]=4*dndmu(me/T); a[4]=4*dndmu(mm/T); a[5]=4*dndmu(ml/T); 
// other solution
    a[0]=a[1]=a[2]=2*np(0);
    a[3]=4*np(me/T); a[4]=4*np(mm/T); a[5]=4*np(ml/T); 


//  Eq.15 :     
    A[0   ]= a[0];   A[3   ]= a[3];                   C[0]=L[0];   // definition of Le
    A[8+1 ]= a[1];   A[8+4 ]= a[4];                   C[1]=L[1];
    A[16+2]= a[2];   A[16+5]= a[5];                   C[2]=L[2]; 
//  Eq. 16    
    A[24  ]= 1;       A[24+3]=-1;                A[24+6]=-1;  // e+nu_e^bar -> d,u^bar 
    A[32+1]= 1;       A[32+4]=-1;                A[32+6]=-1;
    A[40+2]= 1;       A[40+5]=-1;                A[40+6]=-1;

// Eq.18
                                                 A[48+6]=ChiQB(T); A[48+7]=ChiBB(T);  C[7]=dB;
// Eq.19    charge conservation
    A[56+3]=-a[3];   A[56+4]=-a[4];  A[56+5]=-a[5];   A[56+6]=ChiQQ(T); A[56+7]=ChiQB(T); 
/*
double AA[64], CC[8];
for(int i=0;i<8;i++) CC[i]=C[i]; for(int i=0;i<64;i++) AA[i]=A[i]; 
*/
    int err=solveLinEq(8,A, C);

printf(" C[0]=%E   C[2]=%E\n",C[0],C[3]);
printf(" C[1]=%E   C[4]=%E\n",C[1],C[4]);
   
    for(int i=0;i<3;i++) muN[i]=C[i];
    for(int i=3;i<6;i++) muL[i-3]=C[i];
/*
double CCC[8];
      
    for(int i=0;i<8;i++) 
    { CCC[i]=0;
      for(int j=0;j<8;j++) CCC[i]+=AA[i*8+j]*C[j];
      printf("CCC[%d]=%E CC=%E\n", i,CCC[i],CC[i]);
    }  
*/    
    return err; 
}

#endif




static double*epsGrid=NULL;
static double*epsWeight=NULL;

static  double *lnTgrid=NULL;
static  double *fT=NULL;
 int NNT=0;


static double epsDist(double e, void * dat) 
{  
    double c=pow(hEff(0)/hEff(exp(lnTgrid[NNT-1])),1./3.);       
//    c=1;
    e*=c;
    double *d=dat;
    if(e<epsGrid[epsN-2]) return polint1(e,epsN,epsGrid, (double *) dat)*c;
    return d[epsN-2]*exp( -(e-epsGrid[epsN-2])/(epsGrid[epsN-2]-epsGrid[epsN-3])*log(d[epsN-3]/d[epsN-2]))*c;
}


 typedef struct { double T0, // minimal temperature for s-neutrino production 
                        eps, // p/T0  at T=T0
                      sin22, // Theta^2
                        dm2, // ms^2 
                        eta; // (chemical potential)                           ; 
                        int kL, inv; 
                }  nuFIpar;

 nuFIpar AllNuPar;  
 static int n_stat;
 
double LasymmS(double T) 
 {  
    if(fT[0]==0) return 0;
    double res;
    if(n_stat>=0)
    { 
       int l=n_stat; if(l>4) l=3;
//       l=2;
       res=polintN(log(T),l, lnTgrid+n_stat-l+1 , fT+n_stat-l+1);
    } else
    {  double lnT=log(T);
       if(T>lnTgrid[0]) return NAN;
       if(T<lnTgrid[NNT-1]) return NAN;
       res=polint3( log(T), NNT,lnTgrid, fT);
    } 
    return res;    
 }


static  double Cprod(double p,double T, int k, double  L) 
{   

     double GF=1.166E-5; // GeV^{-2}
     double MZ=91.187,MW=80.379; 
     double xl;
     switch(k) 
     { case 1: xl=5.11E-4/T; break;
       case 2: xl=105E-3/T;  break;
       case 3: xl=1.777/T;   break;
     }  

     nuFIpar*par=&AllNuPar;
     double Ll[3]={0,0,0}, dNnu[3],dNch[3], m[3]={5.11E-4,0.105,1.77};
     for(int i=0;i<3;i++) Ll[i]=Lds[i]*(2*M_PI*M_PI/45.*hEff(T)); 
     Ll[k-1]=L;
        
     for(int i=0;i<3;i++) DeltaN_1M[i]=deltaN_1(m[i]/T); 
int lmax=1;
      if(lmax==1) for(int i=0;i<3;i++) DeltaN_3M[i]=0; else for(int i=0;i<3;i++) DeltaN_3M[i]=deltaN_3(m[i]/T);    
     
       // initial values  
// findAsymm2(T, Ll, chiN,chiL);      
     double chi[6]={0,0,0};
     
     for(int l=0;l<lmax;l++)
     { 
     
       double chiOld[6], alpha=0.5;
       for(int i=0;i<6;i++) chiOld[i]=chi[i];
       int err=findAsymm2(T, Ll, chi);
       if(l==0) 
       { double Q=0,dQ=0,dL[3]={0}; 
         for(int i=0;i<3;i++) { Q+=DeltaN_1M[i]*chi[i+3]; dQ+= DeltaN_3M[i]*pow(chi[i+3],3);} 
         for(int i=0;i<3;i++) dL[i]=pow(chi[i],3)*DeltaN_3_0+pow(chi[i+3],3)*DeltaN_3M[i];
         if(fabs(dQ)<1E-2*fabs(Q) && fabs(dL[AllNuPar.kL-1])<1E-2*fabs(Ll[AllNuPar.kL-1]))  break;
//         else printf("Corrections dQ=%E Q=%e kL=%d  dL=%E L=%E\n",dQ,Q, AllNuPar.kL, dL[AllNuPar.kL-1],Ll[AllNuPar.kL-1]);
       } else         
       {  for(int i=0;i<6;i++)  chi[i]=alpha*chi[i]+(1-alpha)*chiOld[i];
          int m;
          for(m=0;m<6;m++) if( fabs(chi[m]-chiOld[m]) > 1E-3*fabs(chi[m])) break;
          if(m==6) { 
//            printf("lastL=%d\n",l); 
          break;}
       }   
     }
     double ml[3]={5.11E-4,0.105,1.777};
     for(int i=0;i<3;i++) 
     { dNnu[i]= chi[i]*(DeltaN_1_0 + chi[i]*chi[i]*DeltaN_3_0); 
       dNch[i]= chi[i+3]*(DeltaN_1M[i] + chi[i+3]*chi[i+3]*DeltaN_3M[i]);
     }       
//printf("k=%d T=%E  L=%E   Ll= %E %E %E dNu= %E %E %E  dNc=%E %E %E\n",k,T,L, Ll[0],Ll[1],Ll[2], dNnu[0],dNnu[1],dNnu[2],dNch[0],dNch[1],dNch[2]);     
     double sinWq=0.23; 
     double Leff=dNnu[k-1]+dNch[k-1];
     for(int i=0;i<3;i++) Leff+= dNnu[i] + 0.5*dNch[i]; //  -0.5 * 6E-10*2*zeta3/(M_PI*M_PI)*hEff(T)/hEff(0); 
//       for(int i=0;i<3;i++) Leff+= dNnu[i] - 0.5*dNch[i];
//printf("Leff=%E\n",Leff);
//printf("err=%d dNnu=%E dNch=%E  Leff=%E\n",err, dNnu[0],dNch[0],   Leff);     
     double gl=  xl<0.5? 7./8. - xl*xl/0.09*(5.820273e-03): g1eff(xl,-1); 
//     gl= g1eff(xl,-1);  printf("gl=%E\n",gl); 
//     gl=0;
     double B=8*sqrt(2)*GF/3* M_PI*M_PI/30*(2*7./8./MZ/MZ+4*gl/MW/MW);
//     B=0;
     double V=B*pow(T,4)*p  -  sqrt(2)*GF*Leff*par->inv*pow(T,3);   // Dolgov eq.277
//     printf("T=%E  B= %E\n",T, B); 
     
//           0101524 Eq.5.13                                0101524 Eq (5.11)     
//   printf(" 137*par->c2a*GF*GF=%.2E(%.2E) 2*sqrt(2)*zeta3/M_PI/M_PI*GF=%e\n",  137*par->c2a*GF*GF,par->B, 2*sqrt(2)*zeta3/M_PI/M_PI*GF);      
     double ga=(k==1)? 3.56 : 2.5;   
     double Gamma=180*zeta3*ga/(7*pow(M_PI,4))*GF*GF*p*pow(T,4);   // 0202122  eq.291
//printf("T=%e p=%E 1+V*2*p/par->dm2=%E\n",T,p,1+V*2*p/par->dm2);     
//     Gamma*=8;
     double res= 0.25*Gamma*par->sin22/(
//     par->sin22+  
     pow(0.5*Gamma*2*p/par->dm2,2)+ 
     pow(1+V*2*p/par->dm2,2));  // left part of 6.6 (0101524) without (fa-fs);    
//     printf("Cprod(T=%E Gamma=%E sin22=%E p=%E res=%E  \n",T,Gamma,par->sin22,p, res);
     return res;
 }




static  void Tderives(double lnT, double * f , double *df)
{    
   double T=exp(lnT);
   nuFIpar*par=&AllNuPar;

   double p=par->eps*T*pow(hEff(T)/hEff(par->T0),1./3.);
//printf("Tderives: par->L=%E\n", par->L);
//printf("p=%e %e %e %e\n",  par->eps, lnT, pow(hEff(T)/hEff(par->T0),1./3.),p); 

//printf("          T=%E Ls(T)=%E par->inv=%d  rest=%E\n", T,Ls(T), par->inv,( 2*M_PI*M_PI/45.*hEff(T)));                                                             
   double res=(1+hEffLnDiff(T)/3)/(Hubble(T)*T)*Cprod(p,T,par->kL,LasymmS(T)*(2*M_PI*M_PI/45.*hEff(T)));
//   *df=-T*res*(1/(exp(p/T-par->eta)+1) -(*f)); 
    *df=-T*res*(1/(exp(p/T)+1) -(*f)) ;
    
//   printf("Cprof=%E   Tderives(%E %E) =%E\n",Cprod(p,T,par->kL, par->Linit*par->inv*( 2*M_PI*M_PI/45.*hEff(T))),  lnT,*f, *df);
} 
 
  
static double *Parr=NULL, *aParr=NULL;

double darkOmegaNu2(double T0, double Tinit, int kL,  double sin22, double Le,double Lm, double Ll,  double msKeV)
{   
  DeltaN_1_0=deltaN_1(0);
  DeltaN_3_0=deltaN_3(0);
  
    Lds[0]=Le; Lds[1]=Lm;Lds[2]=Ll;

//  eps grid: eps=p/T
    epsGrid=realloc(epsGrid,epsN*sizeof(double));
   
    epsWeight=realloc(epsWeight,epsN*sizeof(double));
    double epsLa;
#ifdef XexpX    
    epsLa=epsMax/log(4*epsN);
 #else    
   epsLa=epsMax/log(epsN/2.);
#endif    
//double sum
    for(int i=0;i<epsN;i++)
    { 
//     epsGrid[i]=epsMin*pow(epsMax/epsMin,i/(epsN-1.));
//     epsWeight[i]=pow(epsGrid[i],3)*(log(epsMax)-log(epsMin))/(epsN-1);
#ifdef XepsX
       double z= sqrt( (i+0.5)/epsN); 
       epsGrid[i]=-epsLa*log(1-z);
       epsWeight[i]=0.5*epsLa/(1-z)/z/epsN;
#else        
       epsGrid[i]=epsLa*log(epsN/(epsN-i-0.5));
       epsWeight[i]=epsLa/(epsN-i-0.5);
#endif       
//        sum+=epsWeight[i]*exp(-epsGrid[i])*epsGrid[i]; printf("eps=%E W=%E sum=%E\n",epsGrid[i], epsWeight[i], sum);  
       epsWeight[i]*=epsGrid[i]*epsGrid[i];    
    }  
//    exit(0);
    
    Parr= realloc(Parr, epsN*sizeof(double));
    aParr=realloc(aParr,epsN*sizeof(double));
     
   lnTgrid=realloc(lnTgrid, sizeof(double)); lnTgrid[0]=log(Tinit);
   fT=realloc(fT,sizeof(double)); fT[0]=Lds[kL-1];
    
   AllNuPar.kL=kL;        
   AllNuPar.sin22=sin22;
   AllNuPar.dm2=msKeV*msKeV*1E-12;

   double*Parr_=malloc(epsN*sizeof(double));
   double*aParr_=malloc(epsN*sizeof(double));
   for(int i=0;i<epsN;i++) { Parr[i]=0; aParr[i]=0; }
   double s0_=2*M_PI*M_PI/45*hEff(T0);

   double Tstep=0.95;
//   Tstep=0.99;
   double Tn=Tinit;
   n_stat=0;
int l;   
   for(;;) 
   { 
printf("Tn=%E\n",Tn);   
     n_stat++;
     int n=n_stat; 
     fT=realloc(fT,(n_stat+1)*sizeof(double));
     lnTgrid=realloc(lnTgrid, (n_stat+1)*sizeof(double));
     NNT=n_stat+1;
     for(;;)
     { Tn=Tstep*exp(lnTgrid[n-1]);
//printf("   n=%d  Tn=%E\n",n, Tn);     
       if(Tn<T0) Tn=T0;
       l= n;
       if(l>4)l=4;
       lnTgrid[n]=log(Tn);                                       
       if(l<=2) fT[n]=fT[n-1]; else fT[n]=polintN(lnTgrid[n],l-1, lnTgrid+n-l+1 , fT+n-l+1);
       if( fabs(fT[n]-fT[n-1])> 0.02*fabs(fT[n-1])) { Tstep=pow(Tstep,0.75); } else break;
     }
//printf("ok\n");     
   
//printf("next point T[%d]=%.3E step=1-%.2E   fT(interpol=%d)=%.3E  \n",n, Tn, 1-Tstep, l-1, fT[n]);

     double dY_=0;
     int kEnd=4;      
     for(int k=0;k<=kEnd;k++)
     { 
//printf("k=%d\n",k);
            AllNuPar.inv=1; 
        for(int i=0;i<epsN;i++)  
        {  Parr_[i]=Parr[i]; AllNuPar.eps=epsGrid[i];
           odeint(Parr_+i,1, lnTgrid[n-1], lnTgrid[n], 1E-5, (lnTgrid[n-1]-lnTgrid[n])/4 , Tderives);    
        }
//printf("one odeint\n");        
        AllNuPar.inv=-1;
        for(int i=0;i<epsN;i++)  
        {  aParr_[i]=aParr[i]; AllNuPar.eps=epsGrid[i]; 
           odeint(aParr_+i,1, lnTgrid[n-1], lnTgrid[n], 1E-5, (lnTgrid[n-1]-lnTgrid[n])/4 , Tderives);   
        }

        double drho=0;
        for(int i=0;i<epsN;i++) drho+=((Parr_[i]-aParr_[i]) - (Parr[i]-aParr[i]))*epsWeight[i];  
        drho*=4*M_PI/pow(2*M_PI,3);
        double dY=drho/s0_; 
// if(k) printf("            dY=%E dY_=%e\n",dY,dY_); 
        if(fabs(dY_-dY)<1E-3*fabs(dY)) break; else  dY_=dY; 
        if(k<kEnd)
        {
           fT[n]=fT[n-1]-dY;         
//           if(fT[n]<0) { fT[n]=0; dY_=fT[n-1]-fT[n];}
         }    
     }

//  if(Tn<0.187 && Tn> 0.184) { char mess[20];sprintf(mess, "T[%d]=%.5e",n,Tn); displayPlot(mess,"e",0,5,0,1,"",0,epsDist,Parr);}
     
//      printf("\n");
      for(int i=0;i<epsN;i++) { Parr[i]=Parr_[i]; aParr[i]=aParr_[i];}
      if(Tn==T0) break;
      if( fabs(fT[n]-fT[n-1] ) < 0.01*fabs(fT[n-1]) && Tstep>0.9) Tstep=pow(Tstep,1.5);   
   }

   double rho;
/*   
   rho=0; for(int i=0;i<epsN;i++) rho+=(Parr[i]-aParr[i])*epsWeight[i];  
   rho*=4*M_PI/pow(2*M_PI,3);   
      
   printf("Y_A=%e\n",rho/s0_);
*/   
   rho=0; for(int i=0;i<epsN;i++) rho+=(Parr[i]+aParr[i])*epsWeight[i]; 
   rho*=4*M_PI/pow(2*M_PI,3); 
//   printf("Y=%e\n",rho/s0_);
   
   n_stat=-10;
   free(Parr_); free(aParr_);
   return 2.742E8*rho/s0_*msKeV*1E-6; 
}    


double fNuSterile(double eps, void * arg) 
{  char*txt=(char*) arg;
   char type[20];
   char*ch=strstr(txt,"nus");
    
   if(!ch) return 0; 
   int p; 
   if(1!=sscanf(ch+3,"%d",&p)) p=0;
   if(ch==txt || ch[-1]==' ') return  pow(eps,p)*(epsDist(eps,Parr)+epsDist(eps,aParr));
   if(ch[-1]=='+') return  pow(eps,p)*epsDist(eps,Parr);
   if(ch[-1]=='-') return  pow(eps,p)*epsDist(eps,aParr);
   return 0;
}



int main(void)
{ 
  double inputTab[11][2]=
  { 
    {0.800E-11,13.0E-5},
    {1.104E-11,10.8E-5},
    {1.523E-11,9.60E-5},
    {2.101E-11,8.85E-5},
    {2.899E-11,8.35E-5},
    {4.000E-11,8.00E-5},
    {5.519E-11,7.73E-5},
    {7.615E-11,7.49E-5},
    {10.506E-11,7.24E-5},
    {14.496E-11,7.00E-5},
    {20.000E-11,6.70E-5}
  };
#ifdef QQ
extern int findAsymm(double T, double * L,double * dNnu,double * dNlch);

double LL[3];

double munue(double T)
{ 
   double L[3]={0.035,-0.035,0};
   for(int i=0;i<3;i++) L[i]=LL[i];
   T/=1000;
   for(int i=0;i<3;i++)   L[i]*=2.*M_PI*M_PI/45.*hEff(T); 
   double muN[3],muL[3];
   findAsymm2(T,L,muN,muL);
   return muN[0];  
}

double mue(double T)
{ 
   double L[3]={0.035,-0.035,0};
   for(int i=0;i<3;i++) L[i]=LL[i];
   T/=1000;
   for(int i=0;i<3;i++)   L[i]*=2.*M_PI*M_PI/45.*hEff(T); 
   double muN[3],muL[3];
   findAsymm2(T,L,muN,muL);
   return  muL[0];  
}

double munum(double T)
{ 
   double L[3]={0.035,-0.035,0};
   for(int i=0;i<3;i++) L[i]=LL[i];
   T/=1000;
   for(int i=0;i<3;i++)   L[i]*=2.*M_PI*M_PI/45.*hEff(T); 
   double muN[3],muL[3];
   findAsymm2(T,L,muN,muL);
   return muN[1];  
}

double mum(double T)
{ 
   double L[3]={0.035,-0.035,0};
   for(int i=0;i<3;i++) L[i]=LL[i];
   T/=1000;
   for(int i=0;i<3;i++)   L[i]*=2.*M_PI*M_PI/45.*hEff(T); 
   double muN[3],muL[3];
   findAsymm2(T,L,muN,muL);
   return muL[1];  
}

double munul(double T)
{ 
   double L[3]={0.035,-0.035,0};
   for(int i=0;i<3;i++) L[i]=LL[i];
   T/=1000;
   for(int i=0;i<3;i++)   L[i]*=2.*M_PI*M_PI/45.*hEff(T); 
   double muN[3],muL[3];
   findAsymm2(T,L,muN,muL);
   return muN[2];  
}

double mul(double T)
{ 
   double L[3]={0.035,-0.035,0};
   for(int i=0;i<3;i++) L[i]=LL[i];
   T/=1000;
   for(int i=0;i<3;i++)   L[i]*=2.*M_PI*M_PI/45.*hEff(T); 
   double muN[3],muL[3];
   findAsymm2(T,L,muN,muL);
   return muL[2];  
}



LL[0]=0.035;
LL[1]=-0.035;
LL[2]=0;

double L[3]={-1E-3,1E-3,0}, dNnu[3],dNlch[3];

/*
findAsymm(0.01, L,dNnu, dNlch);
printf("dNu=   %E %E %E\n",dNnu[0],dNnu[1],dNnu[2]);
printf("dNlch= %E %E %E\n",dNlch[0],dNlch[1],dNlch[2]);

exit(0);
*/
/*
double ff( double pi, void * xx)
{ double x= *(double *)xx;
  return pi*pi/(1+exp(pi -x));
}  


double nmu(double x)
{  return simpson_arg(ff, &x,  0, 100, 1E-3,NULL);
} 

displayPlot("n(x)","x",0,4,0,1,"",0,nmu,NULL);
*/

displayPlot("mu/T","T[MeV]",10,5000,1,6,"nue",0,munue,NULL,"e",0,mue,NULL
                                        ,"num",0,munum,NULL,"m",0,mum,NULL
                                        ,"nul",0,munul,NULL,"l",0,mul,NULL );

LL[0]=0.1;
LL[1]=-0.1;
LL[2]=0;

displayPlot("mu/T","T[MeV]",10,5000,1,6,"nue",0,munue,NULL,"e",0,mue,NULL
                                        ,"num",0,munum,NULL,"m",0,mum,NULL
                                        ,"nul",0,munul,NULL,"l",0,mul,NULL );
                                             

exit(0);
#endif

//#define NEW

#ifdef NEW
     double Omega=darkOmegaNu2(0.027, 3, 1, 1.4E-14 ,-0.035, 0.035, 0, 10);

//     double Omega=darkOmegaNu2(0.01, 3, 2, 3E-14 , 0 , 0.018, -0.18, 20);
     double pT=simpson_arg(fNuSterile,"nus 3", 0,10,1E-3,NULL)/simpson_arg(fNuSterile,"nus 2", 0,10,1E-3,NULL);
     printf("Omega h^2=%.2E <p/T>=%.2E\n",  Omega,pT);
     displayPlot("L/s(T)","T",0.01,3,1,1,"",0,LasymmS,NULL);
     displayPlot("e^2f(e)","e=p/T",0,10,0,1,"",0,fNuSterile,"nus 2");
exit(0);
#endif

  for(int n=0;n<1;n++)
  {  double T0=0.01,Tinit=3;  // [GeV]

     double Omega=darkOmegaNu2(T0, Tinit, 2, inputTab[n][0], 0 , inputTab[n][1],0,7.1);
     double pT=simpson_arg(fNuSterile,"nus 3", 0,10,1E-3,NULL)/simpson_arg(fNuSterile,"nus 2", 0,10,1E-3,NULL);
     printf("sin22=%.3E  L/s(init)=%.3E Omega h^2=%.2E <p/T>=%.2E\n",
     inputTab[n][0], inputTab[n][1], Omega,pT);

     displayPlot("L/s(T)","T",0.01,3,1,1,"",0,LasymmS,NULL);
     displayPlot("e^2f(e)","e=p/T",0,10,0,1,"",0,fNuSterile,"nus 2");    
  }

  killPlots();
}
