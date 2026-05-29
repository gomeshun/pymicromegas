/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
#define MASSES_INFO      
  /* Display information about mass spectrum  */
  
#define CONSTRAINTS 
//#define SMODELS
//#define HIGGSBOUNDS 
//#define HIGGSSIGNALS
//#define LILITH

#define OMEGA
  
//#define CDM_NUCLEON     
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

#define CDM_NUCLEUS     
     // Calculate  exclusion rate for direct detection experiments Xenon1T, DarkSide50, CRESST, and PICO 
         
//#define NEUTRINO    
 /*  Neutrino signal of DM annihilation in Sun and Earth */



#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"



int main(int argc,char** argv)
{  int err;

#define NN 100

   if(argc==1) // no file to read. A general analysis 
   { double dM[NN];  
     double DMcm0=1E-17;
     double Mmin=1,Mmax=1000;      
     assignValW("milliQ",0);
     assignValW("DMcm", DMcm0);  
        
     for(int i=0;i<NN;i++) 
     {  
       assignValW("Mdm",Mmin*pow(Mmax/Mmin, (i+0.5)/NN)) ;    
       for(int k=0;k<20;k++)
       {   double Xf, Beps=1E-5; 
           int fast=0;
           assignValW("DMcm",DMcm0);
           sortOddParticles(NULL);
           double Omega=darkOmega(&Xf,fast,Beps,&err);
   
           if(fabs(Omega-0.12) < 0.005) 
           { dM[i]=DMcm0; //printf("M=%E k=%d\n",findValW("Mdm"), k);  
             break;   
           } else  DMcm0*=pow(Omega/0.12, 0.4);    
       }
     }
     displayPlot("d_M for Omegah^2=0.12","Mdm[GeV]",Mmin,Mmax,1,1,"d_M[cm]",NN,dM,NULL);  
     
          
     printf("Direct Detection limit\n");

     Mmin=0.2,Mmax=1000;
     DMcm0=1E-20;
     assignValW("DMcm",DMcm0);
     double CRESST[NN],DS[NN],LZ5T[NN],XENON1T[NN];
     for(int i=0;i<NN;i++){CRESST[i]=NAN; DS[i]=NAN; LZ5T[i]=NAN;}
     int DSc=1,CRESSTc=1,LZ5Tc=1;  // to continue

     for(int i=0;i<NN;i++) 
     { double Mdm=Mmin*pow(Mmax/Mmin,(i+0.5)/NN);   
       assignValW("Mdm",Mdm);
       sortOddParticles(NULL); 
       char *expName;
       double f=DD_factor(AllDDexp, 0.1, Maxwell,&expName);
       double DM=DMcm0*sqrt(f);

       if(strcmp(expName,"DarkSide_2018")==0)  { DS[i]=DM; }       else 
       if(strcmp(expName,"CRESST_2019")==0)    { CRESST[i]=DM; }   else
       if(strcmp(expName,"LZ5T_2024")==0)      {  LZ5T[i]=DM;  }   

       
       if(i) 
       { 
         if(DSc     && isfinite(DS[i-1])     && !isfinite(DS[i])    ) { DS[i]=DM;     DSc=0;    } 
         if(CRESSTc && isfinite(CRESST[i-1]) && !isfinite(CRESST[i])) { CRESST[i]=DM; CRESSTc=0;} 
         if(LZ5Tc   && isfinite(LZ5T[i-1])   && !isfinite(LZ5T[i])  ) { LZ5T[i]=DM;   LZ5Tc=0;  }
       }  
     }    
          
     displayPlot("magnetic dipole M[cm] 90\% exclusion", "Mdm[GeV]",Mmin,Mmax,1,3,"CRESST",NN,CRESST,NULL
                                                                                ,"DS50",  NN,DS,    NULL
                                                                                ,"LZ5T",  NN,LZ5T,  NULL
                                                                       );        
                                                                       
    killPlots();
    return 0;
  }

  char cdmName[20];
  int spin2, charge3,cdim;
                             
  err=readVar(argv[1]);
  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
           
  

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

  printf("milli charge=%.2E  magnetic dipole [cm]=%.2E  Mdm=%.2E\n", findValW("milliQ"), findValW("DMcm"),pMass(CDM[1]));  


  for(int k=1;k<=Ncdm;k++) 
  { 
     qNumbers(CDM[k], &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM[k],  spin2,McdmN[k]); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double csLim;
  if(Zinvisible()) printf("Excluded by Z->invizible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q q-\\bar  Cross Section= %.2E pb\n",csLim);
}
#endif

#ifdef SMODELS
{
  int combineSRs=0;
  char* combineAnas=NULL;//"ATLAS-SUSY-2018-05-ewk,ATLAS-SUSY-2019-08,ATLAS-SUSY-2019-09";
  char*ExcludeBSM=NULL;
  int status=0, smodelsOK=0; 
  double Rvalue, Rexpected, SmoLsig, SmoLmax, SmoLSM;
  double CombRvalue, CombRexpected, CombSmoLsig, CombSmoLmax, CombSmoLSM;

  char analysis[50]={},topology[100]={},smodelsInfo[100];
  char CombAnalyses[200]={};
  int LHCrun=LHC8|LHC13;  //  LHC8 - 8TeV; LHC13 - 13TeV;   
//  int LHCrun=LHC13;  //  LHC13 - 13TeV only;   

  printf("\n\n=====  LHC constraints with SModelS  =====\n\n");

#include "../include/SMODELS.inc" // SLHA interface with SModelS

  printf("SModelS %s \n",smodelsInfo);
  if(smodelsOK) 
  { printf("\n highest r-value = %.2E",Rvalue); 

    if(Rvalue>0) 
    { printf(" from %s, topology: %s ",analysis,topology);
      if(Rexpected>0) 
      { printf("\n expected r = %.2E ",Rexpected);
        if(SmoLsig!=INFINITY)
        { printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.2E %.2E",
                  2*(SmoLsig-SmoLmax), 2*(SmoLSM-SmoLmax) ); 
        }                  
      }
    }  
    if(status==1) { printf("\n excluded by SMS results"); }
    else if(status==0) printf("\n not excluded"); 
    else if(status==-1) printf("\n not not tested by results in SModelS database"); 
    printf("\n");

    // r-value and likelihoods from analysis cvombination
    if(CombRvalue>0) 
    { printf("\n Combination of %s",CombAnalyses);
      printf("\n r-value = %.2E (expected r = %.2E)",CombRvalue, CombRexpected); 
      if(CombRvalue>=1) printf("  --> excluded"); 
      else printf("  --> not excluded"); 
      printf("\n -2log (L_signal/L_max, L_SM/L_max) = %.2E %.2E \n\n", 
                    2*(CombSmoLsig-CombSmoLmax),2*(CombSmoLSM-CombSmoLmax));                     
    }

  } else system("cat smodels.err"); // problem: see smodels.err
}   

#endif 



#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0=3, NHch=1; // number of neutral and charged Higgs particles.
   int HB_id[3]={0,0,0},HB_result[3];
   double  HB_obsratio[3],HS_observ=-1,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50]; 
   NH0=hbBlocksMO("HB.in",&NHch); 
//    NH0= hbBlocksMDL("HB.in",&NHch); 
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("HiggsBounds(%s)\n", HB_version);
   for(int i=0;i<3;i++) if(HB_id[i]) printf("  id= %d  result = %d  obsratio=%.2E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
#endif 
#ifdef HIGGSSIGNALS
   if(HS_observ>=0)
   {
     printf("HiggsSignals(%s)\n",HS_version); 
     printf("  Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_observ,HS_chi2, HS_pval);
   }
#endif   
}
#endif

#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char Lilith_version[50];
   if(LilithMO("Lilith_in.xml"))
   {        
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n",
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }  
   } else printf("LILITH: there is no Higgs candidate\n");
}     
#endif




#ifdef OMEGA
{ printf("\nCalculation of  DM relic density\n");
  
  if(COSMO==FreezeOut) 
  {  printf("    FreezeOut scenario\n");

    int fast=1;
    double Beps=1.E-4, cut=0.01;
    double Omega;  
    int i,err; 
    printf("\n==== Calculation of relic density =====\n");   

    if(Ncdm==1) 
    {  double Xf;
       Omega=darkOmega(&Xf,fast,Beps,&err);
       printf("Xf=%.2e Omega=%.2e err=%d\n",Xf,Omega,err);
       if(err<32)printChannels(Xf,cut,Beps,1,stdout);
    } else
    if(Ncdm==2)
    {
      Omega= darkOmega2(fast,Beps,&err);
      printf("Omega_1h^2=%.2E Omega_2h^2=%.2E err=%d \n", Omega*fracCDM[1], Omega*fracCDM[2],err);
    }else 
    {  
      Omega=darkOmegaN(fast,Beps,&err);
      printf("!Omega=%.2E\n",Omega);
      for(int k=1;k<=Ncdm;k++) printf("   Omega_%d=%.2E\n",k,Omega*fracCDM[k]); 
    }    
  }
  if(COSMO==FreezeIn)  
  { printf("FreezeIn scenario\n");

    double TR=1E6;
    double omegaFi;  
    toFeebleList(CDM[1]);
    VWdecay=0; VZdecay=0;   
    omegaFi=darkOmegaFi(TR,CDM[1],&err);
    printf("omega freeze-in=%.3E\n", omegaFi);
    printChannelsFi(0,0,stdout);
  }
}  
#endif


#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */



//  if(SpA)

#ifdef Dwarf
    if(DwarfSignal(sigmaV,SpA)>1)  printf("Excluded at 95%% confidence level (2002.01229)\n");
#endif   

   if( PlanckCMB(sigmaV,SpA, SpE)>1 ) printf("Excluded at 95%% confidence level (1506.03811)\n"); 

  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }

//  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
//  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  
#endif

#ifdef RESET_FORMFACTORS
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");
 
  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);
//                    To restore default form factors of  version 2  call 
     calcScalarQuarkFF(0.553,18.9,55.,243.5);


  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors  current version  call 
//  calcScalarQuarkFF(0.56,20.2,34,42);


}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        
  double csSIp1,csSIn1,csSDp1,csSDn1, csSIp1_,csSIn1_,csSDp1_,csSDn1_;
  double csSIp2,csSIn2,csSDp2,csSDn2, csSIp2_,csSIn2_,csSDp2_,csSDn2_;
printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

  for(int k=1;k<=Ncdm;k++) 
  {  
    nucleonAmplitudes(CDM[k], pA0,pA5,nA0,nA5);
    printf("%s[%s]-nucleon micrOMEGAs amplitudes\n",CDM[k],antiParticle(CDM[k]));
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*McdmN[k]/(Nmass+ McdmN[k]),2.);
    csSIp1=  SCcoeff*pA0[0]*pA0[0];  csSIp1_=  SCcoeff*pA0[1]*pA0[1];
    csSDp1=3*SCcoeff*pA5[0]*pA5[0];  csSDp1_=3*SCcoeff*pA5[1]*pA5[1];
    csSIn1=  SCcoeff*nA0[0]*nA0[0];  csSIn1_=  SCcoeff*nA0[1]*nA0[1];
    csSDn1=3*SCcoeff*nA5[0]*nA5[0];  csSDn1_=3*SCcoeff*nA5[1]*nA5[1];
                    
    printf("%s[%s]-nucleon cross sections[pb]:\n",CDM[k],antiParticle(CDM[k]));
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIp1,csSIp1_,csSDp1,csSDp1_);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIn1,csSIn1_,csSDn1,csSDn1_);     
  }
}
#endif
  
#ifdef CDM_NUCLEUS
{ char* expName; 
  printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
       if(pval<0.1 )  printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval)); 
  else printf("Not excluded by DD experiments  at 90%% level \n");  
}
#endif 

#ifdef NEUTRINO
if(!CDM[1] || !CDM[2])
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;
  
 printf("\n===============Neutrino Telescope=======  for  "); 
 if(forSun) printf("Sun\n"); else printf("Earth\n");  

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0, 2,"dnu/dE",0,SpectdNdE,nu,"dnu_bar/dE",0,SpectdNdE,nu_bar);
#endif
{ 
    printf(" E>%.1E GeV neutrino flux       %.2E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,nu_bar,NULL));  
} 
  
/* Upward events */
  
  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS
  displayPlot("Upward muons[1/Year/km^2/GeV]","E",Emin,Mcdm/2, 0,1,"mu",0,SpectdNdE,mu);  
#endif
    printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));
  
/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS 
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E",Emin,Mcdm,0,1,"",0,SpectdNdE,mu); 
#endif
  printf(" E>%.1E GeV Contained muon flux %.2E [1/Year/km^3]\n",Emin,spectrInfo(Emin/Mcdm,mu,NULL));
}        
#endif 


#ifdef DECAYS
{ char*  pname = pdg2name(25);
  txtList L;
  double width; 
  if(pname)
  { 
    width=pWidth(pname,&L);  
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  } 
  
}            
#endif

#ifdef CROSS_SECTIONS
{
  char* next,next_;
  double nextM;
    
  next=nextOdd(1,&nextM); 
  if(next && nextM<1000)  
  { 
     double cs, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     char*next_=antiParticle(next);
     Qren=Qfact=nextM; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
     Qren=Qfact;
     cs=hCollider(Pcm,1,nf,Qren, Qfact, next,next_,pTmin,1);
     printf("Production of 'next' odd particle: cs(pp-> %s,%s)=%.2E[pb]\n",next,next_, cs);
  }  
}
 
#endif 

#ifdef CLEAN
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
  system("rm -f   smodels.*");  
#endif 
  
  
  killPlots();  
  return 0;
}
