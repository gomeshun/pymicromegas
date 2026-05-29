#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"



double  h_int( double p, void * arg) 
{  double *argTmmu= arg;
   double T=argTmmu[0], m=argTmmu[1],mu=argTmmu[2];
   double f=1/(exp(sqrt(p*p+m*m)/T - mu/T) -1); 
   return (-f*log(f) +(1-f)*log(1-f))*p*p/M_PI/M_PI/2;
}


int main(void) 
{  double arg[3]={2,0,0};
   double s=simpson_arg(h_int, arg,0,100,1E-3,NULL)/(2*M_PI*M_PI*pow(arg[0],3)/45);
   printf("geff: photons  %E electrons %E \n", g1eff(0,0,1), g1eff(0,0,-1));
   printf("heff: photons  %E electrons %E  \n",h1eff(0,0,1), h1eff(0,0,-1));
}  