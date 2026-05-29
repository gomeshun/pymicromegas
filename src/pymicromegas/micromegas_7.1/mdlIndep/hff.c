#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"

double GF=1.166E-5;
//double MtPole=172.56;

double GammaQpole(double Mh, double *MqPole)
{
   return 3*GF*Mh*pow(*MqPole,2)/(4*sqrt(2)*M_PI)*pow(1- pow(2*(*MqPole)/Mh,2),3./2.);
}

double GammaMtEff(double Mh)
{
   double Mt=MtEff(Mh);  
   return 3*GF*Mh*pow(Mt,2)/(4*sqrt(2)*M_PI)*pow(1- pow(2*Mt/Mh,2),3./2.);
}


double GammaMtbe(double Mh)
{
   double Mt=MtEff(Mh);  
   if(Mt>MtPole) Mt=MtPole;
   double be=pow(1- pow(2*MtPole/Mh,2),2);

   Mt=MtPole*(1-be)+be*Mt;
   return 3*GF*Mh*pow(Mt,2)/(4*sqrt(2)*M_PI)*pow(1- pow(2*Mt/Mh,2),3./2.);
}


double GammaMbEff(double Mh)
{
   double Mb=MbEff(Mh);  
   return 3*GF*Mh*pow(Mb,2)/(4*sqrt(2)*M_PI)*pow(1- pow(2*Mb/Mh,2),3./2.);
}


double GammaMbbe(double Mh)
{
   double Mb=MbEff(Mh);  
   if(Mb>MbPole) Mb=MbPole;
   double be=pow(1- pow(2*MbPole/Mh,2),2);

   Mb=MbPole*(1-be)+be*Mb;
//   printf("Mh=%E MbPole=%E Mb=%E\n", Mh,MbPole,Mb);
   return 3*GF*Mh*pow(Mb,2)/(4*sqrt(2)*M_PI)*pow(1- pow(2*Mb/Mh,2),3./2.);
}




int main(void)
{
/*
      double MtPole=172.56;
static double MbMb=4.183;
static double McMc=1.273;
static double alphaMZ=0.1181;
*/ 
  
  initQCD(0.1181, 1.273, 4.183, MtPole);
  printf("MbPole=%E McPole=%E\n", MbPole, McPole);
  
  printf("MbRun(200)=%E\n",MbRun(200));
  
//  exit(0);

  extern double hWidthCoeff(double Mh, int nl, double *C1);   
  double C1[3];  
  double wh=hWidthCoeff(2.24, 3,C1);
  printf("wh=%E\n", wh);

//exit(0);
  
   displayPlot("Gamma h->t,T","Mh",2*MtPole,500,0,3,"PoleMass",0,GammaQpole,&MtPole
                                                    ,"MtEff", 0,GammaMtEff, NULL
                                                    ,"Mtbe",  0,GammaMtbe, NULL );
 
  
  displayPlot("Gamma h->t,T","Mh",500,10000,1,3,"PoleMass",0,GammaQpole,&MtPole
                                                    ,"MtEff",   0,GammaMtEff, NULL 
                                                    ,"Mtbe",    0,GammaMtbe, NULL);
  


  displayPlot("Gamma h->b,B","Mh",2*MbPole,20,0,3,"PoleMass", 0,GammaQpole,&MbPole
                                                    ,"MbEff", 0,GammaMbEff, NULL
                                                    ,"Mbbe",  0,GammaMbbe,   NULL );

  displayPlot("Gamma h->b,B","Mh",20,200,1,3,"PoleMass", 0,GammaQpole,&MbPole
                                                    ,"MbEff", 0,GammaMbEff, NULL
                                                    ,"Mbbe",  0,GammaMbbe,   NULL );
                                                                                                         
                                                    
}