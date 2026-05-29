#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

int err;

typedef   struct{ double mu, chi; int  eta;}    geff_struct;


static double g1eff_int(double x,void *par_)
{  if(x==0 || x==1) return 0;
   geff_struct*par=par_;
   double p=-(3+par->mu/4)*log(x);
   double dpdx=(3+par->mu/4)/x;
   double e=sqrt(p*p+par->mu*par->mu)-par->chi;
   double den;

//   if(par->eta==-1 && e<1E-8) den=e*(1+e/2); else 
   den=exp(e)+par->eta;

   return  dpdx*e*p*p/den/M_PI/M_PI;
}


static double h1eff_int(double x,void *par_)
{  if(x==0 || x==1) return 0;
   geff_struct*par=par_;
   double p=-(3+par->mu/4)*log(x);
   double dpdx= (3+par->mu/4)/x;
   double e=sqrt(p*p+par->mu*par->mu)-par->chi;
   double den;

   den=exp(e)+par->eta;

   double f=1/den;
   if(par->eta<0)  return -dpdx*( f*log(f) - (1+f)*log(1+f))*p*p/(2*M_PI*M_PI) ; // bosons
   
    else           return  dpdx*(-f*log(f) - (1-f)*log(1-f))*p*p/(2*M_PI*M_PI) ; // fermions
//   return  (par->mu*par->mu+4./3.*p*p)*p*p/den/x/M_PI/M_PI/e;
}


double g1eff_(double mu, double chi,int eta)
{ if(!isfinite(mu)) return 0;
  geff_struct par;
  par.eta=-eta;  //!!!
  par.mu=mu;
  par.chi=chi;
  int err;
  double r= simpson_arg(g1eff_int,&par, 0, 1,1E-4,NULL)/(M_PI*M_PI/30);
  return r;  
}


double h1eff_(double mu,double chi, int eta)
{ if(!isfinite(mu)) return 0;
  geff_struct par;
  par.eta=-eta;  //!!!
  par.mu=mu;
  par.chi=chi;
  int err;
  double r= simpson_arg(h1eff_int,&par, 0, 1,1E-4,&err)/(2*M_PI*M_PI/45);
  if(err) 
  {
     displayPlot("h1eff_int","x",0,1,0,1, "",0,h1eff_int,&par);
     exit(0);
  
  }
  
  return r;
}



static double DNintegrand(double x, void* arg)
{ if(x==0) return 0;
  double* arg_=(double*) arg;
  double eps= arg_[0], chi=arg_[1];
  double p=-2*log(x);
//  dp=2 dx/x 
  return 2*p*p/x*(1/(1+exp(sqrt(eps*eps+p*p)-chi)) - 1/(1+exp(sqrt(eps*eps+p*p)+chi)));   
} 

static double Nintegrand(double x, void* arg)
{ if(x==0) return 0;
  double* arg_=(double*) arg;
  double eps= arg_[0], chi=arg_[1];
  double p=-2*log(x);
//  dp=2 dx/x 
  return 2*p*p/x*(1/(1+exp(sqrt(eps*eps+p*p)-chi)) + 1/(1+exp(sqrt(eps*eps+p*p)+chi)));   
} 

double Ntot(double eps, double chi)
{ double v[2];
  v[0]=eps; v[1]=chi; 
  return 1/(2*M_PI*M_PI)*simpson_arg(Nintegrand, v,0,1,1E-5,NULL);
}



double DeltaN(double eps, double chi)
{ double v[2];
  v[0]=eps; v[1]=chi; 
  return 1/(2*M_PI*M_PI)*simpson_arg(DNintegrand, v,0,1,1E-5,NULL);
}

double DeltaN_(double chi, double* eps) { return DeltaN(*eps,chi);} 

double DeltaN_approx(double chi, double*arg)  { return  chi*(arg[0]+ chi*chi*arg[1]); }

static  double phi(double x) { double c= 2*cosh(sqrt(x)/2); return 1/(c*c);  } 

/*
static double DN1integrand(double x, void* arg)
{ if(x==0) return 0;
  double* arg_=(double*) arg;
  double eps= arg_[0];
  double p=-(1+eps/4)*log(x);
//  dp=2 dx/x 
  return (1+eps/4)*p*p/x*phi(eps*eps+p*p) ;
}
*/

static double DN1integrand(double x, void* arg)
{ if(x==0) return 0;
  double* arg_=(double*) arg;
  double eps= arg_[0];
  double p= -(3+eps/4)*log(x);
  double dpdx=(3+eps/4)/x; 
  return  dpdx*p*p*phi(eps*eps+p*p) ;
}


double DeltaN_1(double eps)  { 
//if(eps>10) return 0;

 double r= simpson_arg(DN1integrand, &eps,0,1,1E-4,&err)/(M_PI*M_PI);


//if(err) 
/*
{  displayPlot("DeltaN_1","x",0,1,0,1,"",0,DN1integrand, &eps);
   printf("eps=%E\n",eps); 
   exit(0);
}
*/
   return r;

}

static double DN3integrand(double x, void* arg)
{ if(x==0) return 0;
  double* arg_=(double*) arg;
  double eps= arg_[0];
  double p=-2*log(x);
//  dp=2 dx/x 
  return 2*p*p/x*phi(eps*eps+p*p)*(1-6*phi(eps*eps+p*p)) ;
}

double DeltaN_3(double eps)  { //if(eps>10) return 0;
return simpson_arg(DN3integrand, &eps,0,1,1E-4,&err)/(6*M_PI*M_PI);}
 

double dNdivN(double eps) { return  DeltaN(eps,2)/Ntot(eps,2);}  

#define zeta3 1.202

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


/*static*/ int findAsymm(double T, double * L, double * chiN,double * chiL)  // T^3 normalized asymmetries   
{   double me=5.11E-4, mm=0.105, ml=1.777;                                   // masses of charged  leptons
    double nBng=6E-10,   dB= nBng*2*zeta3/M_PI/M_PI*hEff(T)/hEff(0);         // barion asymmetry normalized on T^3   

    double A[64]={0},C[8]={0};
//  A[i*8+ j]
//  j=0,1,2 - neutrino asymmetry 
//  j=3,4,5 - charger lepton asymmetry
//  j=6 - quark asymmetry  u-d  
//  j=7 - quark asymmetry 2u+d
    double a[6]; 
    a[0]=a[1]=a[2]=DeltaN_1(0);   // 2 because of  asymmetry.Particles increas, aParticle decreas.
    a[3]=2*DeltaN_1(me/T); a[4]=2*DeltaN_1(mm/T); a[5]=2*DeltaN_1(ml/T); 

//  Eq.15 :     
    A[0   ]= a[0];   A[3   ]= a[3];                   C[0]=L[0]-DeltaN_3(0)*pow(chiN[0],3) - 2*DeltaN_3(me/T)*pow(chiL[0],3) ;   // definition of Le
    A[8+1 ]= a[1];   A[8+4 ]= a[4];                   C[1]=L[1]-DeltaN_3(0)*pow(chiN[1],3) - 2*DeltaN_3(mm/T)*pow(chiL[1],3);
    A[16+2]= a[2];   A[16+5]= a[5];                   C[2]=L[2]-DeltaN_3(0)*pow(chiN[2],3) - 2*DeltaN_3(ml/T)*pow(chiL[2],3);
//  Eq. 16    
    A[24  ]= 1;       A[24+3]=-1;                A[24+6]=-1;  // e+nu_e^bar -> d,u^bar 
    A[32+1]= 1;       A[32+4]=-1;                A[32+6]=-1;
    A[40+2]= 1;       A[40+5]=-1;                A[40+6]=-1;

// Eq.18
                                                 A[48+6]=ChiQB(T); A[48+7]=ChiBB(T);       C[6]=dB;
// Eq.19    charge conservation
    A[56+3]=-a[3];   A[56+4]=-a[4];  A[56+5]=-a[5];   A[56+6]=ChiQQ(T); A[56+7]=ChiQB(T); C[7]=DeltaN_3(me/T)*pow(chiL[0],3)
                                                                                              +DeltaN_3(mm/T)*pow(chiL[1],3)
                                                                                              +DeltaN_3(ml/T)*pow(chiL[2],3);
                                                               //                           C[7]=0;    
    int err=solveLinEq(8,A, C);
//    printf("err=%d\n",err);

    for(int i=0;i<3;i++) chiN[i]=C[i];
    for(int i=3;i<6;i++) chiL[i-3]=C[i];

    return err; 
}


double heL(double T) 
{  double ml[3]={5.11E-4,0.105,1.777}; 
   double h=0;
   double s=0;
   for(int i=0;i<3;i++) s+=4*h1eff_(ml[i]/T,0,-1)+2*h1eff_(0,0,-1); 
   return s;
}   





int main(void)
{
  double Ls[3]={0.035,-0.035,0};

  displayPlot("hEff and hL","T[GeV]", 1E-3, 1E3,1,2,"hEff",0, hEff,NULL, "hL",0,heL,NULL);



#ifdef oneP  //  testing of single point

  double T=3; // [GeV]

  double chiN[3]={0,0,0}, chiL[3]={0,0,0};
  
  for(int i=0;i<3;i++)   Ls[i]*=2.*M_PI*M_PI/45.*hEff(T); 
  for(int k=0;k<1000;k++) 
  { double chiNold[3], chiLold[3];
    for(int i=0;i<3;i++) { chiNold[i]=chiN[i]; chiLold[i]=chiL[i];} 
    findAsymm(T,Ls,chiN,chiL);
    if(k) for(int i=0;i<3;i++) { chiN[i]=0.5*chiN[i]+0.5*chiNold[i]; chiL[i]=0.5*chiL[i]+0.5*chiLold[i];} 
    printf("%E %E %E %E %E %E\n", chiN[0],chiN[1],chiN[2], chiL[0],chiL[1],chiL[2]);
  }
     
exit(0);
#endif

#define N 100    
  double T1=1, T2=1000; // MeV
  double chiNarr[3][N],chiLarr[3][N];
  double hEffArr[N],gEffArr[N]; 

for(int k=0;k<N;k++) 
{ double T=T1*pow(T2/T1, (0.5+k)/N)/1000;
   hEffArr[k]=hEff(T);
}

  
  double alpha=0.1;
  for(int k=0;k<N;k++) 
//  for(int k=0;k<1;k++)
  {  double T=T1*pow(T2/T1, (0.5+k)/N)/1000; // GeV
     hEffArr[k]=hEff(T);
  
     double LT[3];
     double chiN[3]={0,0,0}, chiL[3]={0,0,0};

for(int l=0;l<5;l++)
{    double ml[3]={5.11E-4,0.105,1.777}; 
     for(int i=0;i<3;i++)   LT[i]=Ls[i]*2.*M_PI*M_PI/45.* hEffArr[k];  // normalized on T^3

     for(int l=0;l<300;l++) 
     {  double chiNold[3], chiLold[3];
        for(int i=0;i<3;i++) { chiNold[i]=chiN[i]; chiLold[i]=chiL[i];} 
        findAsymm(T,LT,chiN,chiL);
        if(l) for(int i=0;i<3;i++) { chiN[i]=alpha*chiN[i]+(1-alpha)*chiNold[i]; chiL[i]=alpha*chiL[i]+(1-alpha)*chiLold[i];} 
     } 

     hEffArr[k]=hEff(T);
  
     for(int i=0;i<3;i++) hEffArr[k]-= 2*h1eff_(0,0,-1)-h1eff_(0,chiN[i],-1) -h1eff_(0,-chiN[i],-1) ;     
     for(int i=0;i<3;i++) hEffArr[k]-= 4*h1eff_(ml[i]/T,0,-1) - 2*h1eff_(ml[i]/T,chiL[i],-1) -2*h1eff_(ml[i]/T,-chiL[i],-1);

//for(int i=0;i<3;i++) printf( "chiN=%E h1eff_(0,0,-1) %E => %E %E\n",chiN[i],2*h1eff_(0,0,-1),h1eff_(0,chiN[i],-1), h1eff_(0,-chiN[i],-1));
//for(int i=0;i<3;i++) printf( "chiL=%E h1eff_(ml[i]/T,0,-1) %E => %E %E\n",chiN[i],2*h1eff_(ml[i]/T,0,-1),h1eff_(ml[i]/T,chiN[i],-1), h1eff_(ml[i]/T,-chiN[i],-1));       

}     
     for(int i=0;i<3;i++) { chiNarr[i][k]=chiN[i]; chiLarr[i][k]=chiL[i]; }
     
     printf("T=%E  heff=%E => %E\n", T, hEff(T), hEffArr[k]);

 double ml[3]={5.11E-4,0.105,1.777}; 
     
     gEffArr[k]=gEff(T);
     for(int i=0;i<3;i++) gEffArr[k]-= 2*g1eff_(0,0,-1)-g1eff_(0,chiNarr[i][k],-1) -g1eff_(0,-chiNarr[i][k],-1) ;
     for(int i=0;i<3;i++) gEffArr[k]-= 4*g1eff_(ml[i]/T,0,-1) - 2*g1eff_(ml[i]/T,chiLarr[i][k],-1) -2*g1eff_(ml[i]/T,-chiLarr[i][k],-1);
   
  }     
  
//exit(0);  
printf("QQ\n");
  displayPlot("chi","T[Mev]",T1,T2,1,6,"ne",N,chiNarr[0],NULL
                                      ," e",N,chiLarr[0],NULL
                                      ,"nm",N,chiNarr[1],NULL
                                      ," m",N,chiLarr[1],NULL
                                      ,"nl",N,chiNarr[2],NULL
                                      ," l",N,chiLarr[2],NULL
                                      ); 
                                      
double he(double T) { return hEff(T/1000);}
double ge(double T) { return gEff(T/1000);}

  displayPlot("hEff","T[Mev]",T1,T2,1,2,
       "hEff_orig",0,he,NULL, 
       "hEff_improved",N,hEffArr,NULL);           
    displayPlot("gEff","T[Mev]",T1,T2,1,2,
       "gEff_orig",0,ge,NULL, 
       "gEff_improved",N,gEffArr,NULL);

printf("h1eff_(0,0,-1)=%E g1eff_(0,0,-1)=%E\n", h1eff_(0,0,-1), g1eff_(0,0,-1));

printf("h1eff_(0,0,1)=%E g1eff_(0,0,1)=%E\n", h1eff_(0,0,1), g1eff_(0,0,1));

double ml[3]={5.11E-4,0.105,1.777};                            

  for(int i=0;i<10;i++) printf("%.4E, ",log(DeltaN_3(i)));
  printf("\n"); 
printf("%.4E, %.4E\n", log(DeltaN_3(20)), log(DeltaN_3(40)));  
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

/*  
  displayPlot("delta_1","x",0,40,0,2,"integration",0,DeltaN_1,NULL, "intepolation", 0, deltaN_1,NULL);

  displayPlot("delta_3","x",0,40,0,2,"integration",0,DeltaN_3,NULL, "intepolation", 0, deltaN_3,NULL);
*/
  
}

