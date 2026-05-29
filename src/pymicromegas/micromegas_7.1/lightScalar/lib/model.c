
#include"../../include/micromegas.h"
#include"../../include/micromegas_aux.h"
#include"pmodel.h"

static int getdata(char * f, double ** xx, double ** yy)
{
  FILE*F;
  char fname[200];
  sprintf(fname,"data/%s",f); // for micromegas
  F=fopen(fname,"r");
  if(!F)
  {  sprintf(fname,"../../data/%s",f);  // for CalcHEP
     F=fopen(fname,"r");
  }
  if(!F) { printf(" Can not fild file %s\n",f); return -1;}
  char buff[100];
  fscanf(F,"%[^\n]\n",buff);
  int n=0;
  for(;;)
  {  double x,y;
     if( fscanf(F,"%lf %lf",&x,&y)!=2 ) break;
     n++;
     *xx=realloc(*xx,n*sizeof(double));
     *yy=realloc(*yy,n*sizeof(double));
     (*xx)[n-1]=x;
     (*yy)[n-1]=y;
  }
  fclose(F);
  if(n<=3) return -1;
  return n;
}


static int getdata4(char * f, double ** xx, double ** yy, double ** zz, double ** ww)
{
  FILE*F;
  char fname[200];
  sprintf(fname,"data/%s",f); // for micromegas
  F=fopen(fname,"r");
  if(!F)
  {  sprintf(fname,"../../data/%s",f);  // for CalcHEP
     F=fopen(fname,"r");
  }
  if(!F) { printf(" Can not find file %s\n",f); return -1;}
  char buff[100];
  fscanf(F,"%[^\n]\n",buff);
  int n=0;
  for(;;)
  {  double x,y,z,w;
     if( fscanf(F,"%lf %lf %lf %lf ",&x,&y,&z,&w)!=4 ) break;
     n++;
     *xx=realloc(*xx,n*sizeof(double));
     *yy=realloc(*yy,n*sizeof(double));
     *zz=realloc(*zz,n*sizeof(double));
     *ww=realloc(*ww,n*sizeof(double));
     (*xx)[n-1]=x;
     (*yy)[n-1]=y;
     (*ww)[n-1]=w;
     (*zz)[n-1]=z;
  }
  fclose(F);
  if(n<=3) return -1;
  return n;
}

 double abs_delta_k(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("abs_delta_k.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double abs_gamma_pi(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("abs_gamma_pi.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double arg_delta_k(double m)
  { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("arg_delta_k.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double arg_gamma_pi(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("arg_gamma_pi.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double abs_delta_pi(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("abs_delta_pi.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double abs_theta_k(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("abs_theta_k.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double arg_delta_pi(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("arg_delta_pi.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double arg_theta_k(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("arg_theta_k.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double abs_gamma_k(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("abs_gamma_k.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double abs_theta_pi(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("abs_theta_pi.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double arg_gamma_k(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("arg_gamma_k.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double arg_theta_pi(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double n=0;

   if(n==0) n=getdata("arg_theta_pi.dat",&x,&y);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double G_btpz_pi(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double *z=NULL;
   static double *w=NULL;
   static double n=0;

   if(n==0) n=getdata4("G_btpz_pi.dat",&x,&y,&z,&w);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double G_dgl_pi(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double *z=NULL;
   static double *w=NULL;
   static double n=0;

   if(n==0) n=getdata4("G_dgl_pi.dat",&x,&y,&z,&w);
   if(n<0) return 0;
   printf("m=%.3e gpi=%.3e\n",m,n);
   return polint3(m,n,x,y);
 }

 double G_btpz_k(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double *z=NULL;
   static double *w=NULL;
   static double n=0;

   if(n==0) n=getdata4("G_btpz_k.dat",&x,&y,&z,&w);
   if(n<0) return 0;
   return polint3(m,n,x,y);
 }

 double G_dgl_k(double m)
 { static double *x=NULL;
   static double *y=NULL;
   static double *z=NULL;
   static double *w=NULL;
   static double n=0;

   if(n==0) n=getdata4("G_dgl_k.dat",&x,&y,&z,&w);
   if(n<0) return 0;
   printf("m=%.3e gk=%.3e\n",m,n);
   return polint3(m,n,x,y);
 }
