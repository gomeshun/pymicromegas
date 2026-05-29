#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

   char * name[5]={"micrO", "Pythia_8", "PPPC", "CosmiXs", "CosmiXs2"};
   char * out[7]={"A","e","p","ne","nm","ml","D"};  

double energyTest(double M, char*pdgTxt)
{  double Etot=0;
   double sp[NZ];
   int pdg; 
   sscanf(pdgTxt,"%d",&pdg);   
   for(int Nout=0;Nout<7;Nout++) 
   { basicSpectra(M,pdg,Nout,sp);
     int err=0;
     double xMin=1E-6;
     for(int i=NZ-1;i;i--) if(sp[i]==0) xMin=exp(Zi(i)); else break;
     double xMax=1;
      for(int i=1;i<NZ;i++) if(sp[i]==0) xMax=exp(Zi(i)); else break;          
     double e=simpson_arg((funcArg)eSpectdNdE,sp,M*xMin,M*xMax,1E-3,&err);
     if(Nout==2) e+=0.939*simpson_arg((funcArg)SpectdNdE,sp,  M*xMin,M*xMax,1E-3,&err);
     if(Nout==6) e+=2*0.939*simpson_arg((funcArg)SpectdNdE,sp,M*xMin,M*xMax,1E-3,&err);
/*      
     if(err)
     {  
        char txt[200]; 
        sprintf(txt,"Problem in integration err=%d: SpectraFlag=%s M=%.2E pdg=%d Nout=%s",err, name[SpectraFlag],M, pdg,out[Nout]);
        printf("%s\n",txt);
//        displayPlot(txt,"E",M*1E-6,M,1,1,"EdNdE",0,eSpectdNdE, sp);  exit(1);
     }
*/         
    if(Nout) Etot+=2*e; else Etot=e; 
//    printf("Nout=%d e=%E\n",Nout,e);
   }
   return Etot/(2*M)-1;
}                                         


int main(int argc,char** argv)
{ 
   double M=100;       //  DM mass
   int pdg=11,         //  b-quark  
   Nout=0;             //  photons
   double xMin=1E-5;   //  xMin for plots 
 
   SpectraFlag=0;
   double sp[7][NZ];
   char title[100];
   sortOddParticles(NULL);


   for(SpectraFlag=0;SpectraFlag<4;SpectraFlag++) 
   { basicSpectra(M,pdg,Nout,sp[SpectraFlag]);
     printf("%E\n",energyTest(M,"5"));
   }


   sprintf(title,"%s,%s -> %s",pdg2name(pdg),pdg2name(-pdg),out[Nout]);
   displayPlot(title,"E",xMin*M,M,1,4,name[0],0,eSpectdNdE,sp[0]
                                      ,name[1],0,eSpectdNdE,sp[1]
                                      ,name[2],0,eSpectdNdE,sp[2]
                                      ,name[3],0,eSpectdNdE,sp[3]);
   exit(0);                                   

   int lnx=(xMin<1E-2)? 1:0; 
   lnx=1;


for(Nout=0;Nout<6;Nout++)
{   
   for(SpectraFlag=0;SpectraFlag<5;SpectraFlag++) basicSpectra(M,pdg,Nout,sp[SpectraFlag]); 
   sprintf(title,"EdNdE for  %d,%d -> %s. M=%.3E",pdg,-pdg,out[Nout],M);  
   displayPlot(title, "E", xMin*M,M,lnx,1 
                                       ,name[0],0,SpectdNdE,sp[0] 
                                       ,name[1],0,SpectdNdE,sp[1]
                                       ,name[2],0,SpectdNdE,sp[2]
                                       ,name[3],0,SpectdNdE,sp[3]
                                       ,name[4],0,SpectdNdE,sp[4]
                                       ); 
}                                           
  SpectraFlag=0;
//exit(0);
// Unsectainty  
   spectUncert=1;
   spectraUncertainty(M, pdg, Nout, sp[4]);
   basicSpectra(M,pdg,Nout,sp[5]);
   sprintf(title,"Uncertainty  EdNdE for  %d[A,%d -> %s. M=%.3E",pdg,-pdg,out[Nout],M);    
   displayPlot(title, "E", xMin*M,M,1,2,"delta(EdNdE)",0,eSpectdNdE,sp[4],"EdNdE",0,eSpectdNdE,sp[5]); 
   spectUncert=0;
   
// Test energy conservation                                                      
   SpectraFlag=3;
   
   sprintf(title,"Energy conservation E/(2Mdm)-1 for %s spectra",name[SpectraFlag]);

   double Mmin=SpectraFlag ? 5: 0.6;

  
   displayPlot(title, "Mdm", Mmin  ,  1000  ,1,  4
//                                ,"pi0",0,energyTest,"111"
//                                ,"e",  0,energyTest,"11"
//                                ,"tay",0,energyTest,"15"
//                                ,"W",0,energyTest,"24"
//                                ,"h",0,energyTest,"25"
                                ,"g",0,energyTest,"21"
                                ,"s",0,energyTest,"3"
                                ,"u",0,energyTest,"2"
                                ,"d",0,energyTest,"1"
//                                ,"c",0,energyTest,"4"   
      );

}                                           

