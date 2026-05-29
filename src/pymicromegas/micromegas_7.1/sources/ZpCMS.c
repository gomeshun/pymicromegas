#include<math.h>
#include"micromegas.h"
#include"micromegas_aux.h"


#define Ndat 35

static   double MZp[Ndat]=   {200,           400,          600,          800,          1000,         1200,         1400,         1600,         1800,         2000,         2200,         2400,         2600,         2800,         3000,         3200,         3400,         3600,         3800,         4000,         4200,         4400,         4600,         4800,         5000,         5200,         5400,         5600,         5800,         6000,         6200,         6400,         6600,         6800,         7000        };
static   double cdMax[Ndat]= { 5.133874E-07, 9.847946E-07, 4.171765E-06, 2.968828E-06, 4.771844E-06, 7.573147E-06, 1.198427E-05, 1.287125E-05, 1.778244E-05, 3.428707E-05, 5.033839E-05, 4.938690E-05, 1.188838E-04, 1.689959E-04, 1.582952E-04, 2.552646E-04, 4.711327E-04, 5.389394E-04, 5.776558E-04, 8.204217E-04, 1.238300E-03, 1.956533E-03, 3.034897E-03, 4.831026E-03, 7.672444E-03, 1.232210E-02, 1.987964E-02, 3.038081E-02, 4.464066E-02, 6.575673E-02, 9.671937E-02, 1.419744E-01, 2.091076E-01, 3.082024E-01, 4.534962E-01};
static   double cuMax[Ndat]= { 3.936680E-07, 7.021297E-07, 2.836800E-06, 1.976658E-06, 3.096489E-06, 4.846399E-06, 7.546297E-06, 7.976197E-06, 1.080608E-05, 2.040595E-05, 2.930582E-05, 2.800567E-05, 6.550551E-05, 9.050536E-05, 8.200518E-05, 1.280050E-04, 2.280048E-04, 2.530047E-04, 2.610045E-04, 3.590044E-04, 5.230042E-04, 7.960041E-04, 1.190004E-03, 1.840004E-03, 2.830004E-03, 4.410004E-03, 6.910003E-03, 1.020000E-02, 1.460000E-02, 2.080000E-02, 2.970000E-02, 4.240000E-02, 6.060000E-02, 8.650000E-02, 1.240000E-01};

double fcdMax(double m) { return polint3(m,Ndat,MZp,cdMax);}
double fcuMax(double m) { return polint3(m,Ndat,MZp,cuMax);}


double ZpLimCMS(char *Zp)  // 2103.02708,   1010.6058
{
  int pdg,spin2,charge3,cdim;
  pdg=qNumbers(Zp, &spin2, &charge3, &cdim);
  if(pdg==0 || spin2!=2 || charge3!=0|| cdim!=1) { printf("%s is not a Z'-like particle\n"); return NAN;}
  
  double MP=pMass(Zp);

//printf("%e <    MP=%e  < %e \n",MZp[0],MP,MZp[Ndat-1]);
  
  if(MP<MZp[0] || MP>MZp[Ndat-1]) return 0; // no data

  
  char *d=pdg2name(1),*D=pdg2name(-1),*u=pdg2name(2),*U=pdg2name(-2),*e=pdg2name(11),*E=pdg2name(-11),*m=pdg2name(13),*M=pdg2name(-13);
  if(!d || !D || !u ||!U || !e || !E|| !m || !M) return NAN;

//  printf(" particles: %s %s %s %s %s %s %s %s\n",d,D,u,U,e,E,m,M);

  double cd=0,cu=0;
  
  lVert* ZpQq=getVertex(D,d,Zp,NULL);
  if(ZpQq)
  { 
     double coeff[10];
     getNumCoeff(ZpQq,coeff);
     if(ZpQq->nTerms > 2) { printf("unexpected type of interaction\n"); return NAN; }
     for(int k=0;k<ZpQq->nTerms;k++) 
     if( strcmp(ZpQq->SymbVert[k],"G(m3)")==0)    cd+=2*pow(coeff[k],2); else 
     if( strcmp(ZpQq->SymbVert[k],"G5*G(m3)")==0) cd+=2*pow(coeff[k],2); else 
     { printf("unexpected type of interaction\n"); return NAN; }
  }
  ZpQq=getVertex(U,u,Zp,NULL);
  if(ZpQq)
  { 
     double coeff[10];
     getNumCoeff(ZpQq,coeff);
     if(ZpQq->nTerms > 2) { printf("unexpected type of interaction\n"); return NAN; }
     for(int k=0;k<ZpQq->nTerms;k++) 
     if( strcmp(ZpQq->SymbVert[k],"G(m3)")==0)    cu+=2*pow(coeff[k],2); else 
     if( strcmp(ZpQq->SymbVert[k],"G5*G(m3)")==0) cu+=2*pow(coeff[k],2); else 
     { printf("unexpected type of interaction\n"); return NAN;} 
  }
  txtList L;
  pWidth(Zp,&L);
  char dielectron[20], dimuon[20];
  sprintf(dielectron,"%s,%s",e,E);
  sprintf(dimuon,"%s,%s",m,M);
  double Br=findBr(L,dielectron)+findBr(L,dimuon);
  Br/=2;
//  printf("l-branching: %e\n",Br);    
  cd*=Br;
  cu*=Br;
  double cdM=polint3(MP,Ndat,MZp,cdMax);
  double cuM=polint3(MP,Ndat,MZp,cuMax); 
//  printf("cu=%E cuMax=%E   cd=%e cdMax=%E\n",cu,cuM,cd,cdM);
//  return  erf( 1.386*(cd/cdM+cu/cuM));
    return cd/cdM+cu/cuM;
}   

