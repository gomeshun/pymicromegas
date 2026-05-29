/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
#define MASSES_INFO (1UL << 0)
  /* Display information about mass spectrum  */
  
#define CONSTRAINTS (1UL << 1)

#define MONOJET      (1UL << 2)
#define HIGGSBOUNDS  (1UL << 3)
#define HIGGSSIGNALS (1UL << 4)
#define LILITH       (1UL << 5)
#define SMODELS      (1UL << 6)
  
#define OMEGA    (1UL << 7) /*  Calculate Freeze out relic density and display contribution of  individual channels */
#define FREEZEIN (1UL << 8) /*  Calculate relic density in Freeze-in scenario  */
   
#define INDIRECT_DETECTION  (1UL << 9)
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line 
     of sight; 
     Calculate galactic propagation of positrons and antiprotons.      
  */
      
#define RESET_FORMFACTORS (1UL << 10)
  /* Modify default nucleus form factors, 
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */     
#define CDM_NUCLEON (1UL << 11)    
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

#define CDM_NUCLEUS (1UL << 12)
  /* Calculate number of events for 1kg*day and recoil energy distibution
      for various nuclei
  */
#define NEUTRINO (1UL << 13)   
 /*  Neutrino signal of DM annihilation in Sun and Earth */

#define DECAYS (1UL << 14)

#define CROSS_SECTIONS (1UL << 15)
  
/*===== end of Modules  ======*/

/*===== Options ========*/
#define SHOWPLOTS (1UL << 16)
     /* Display  graphical plots on the screen */ 
#define CLEAN (1UL << 17)
/*===== End of DEFINE  settings ===== */

#include<string.h>

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


int readVar_str(int nVar, char **varnames, char **strvals){
  int err = 0;
  char *errstr;
  double val;
  for(int i=0; i<nVar;i++){
    val = strtod(strvals[i],&errstr);
    if(*errstr != '\0'){
      printf("invalid input:%s in %dth argument\n", errstr,i+1);
      exit(1);
    }
  err=assignVal(varnames[i],val);
  if(err==1){ err=i+1; break; }
  }
  return err;
}



int main(int argc,char** argv)
{  int err = 0;
   char *errstr;
   char cdmName[10];
   int spin2, charge3,cdim;
   unsigned long flags = 0;
   int n_inputvals = 0;
   const int n_options = 3;  // <integer for flags>, <the number of parameterss>, <DOFfilename or None>

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */

  VZdecay=0; VWdecay=0;  

  if(argc<1+n_options)
  { 
      printf("usage:  ./main <integer for flags> <the number of parameters> <DOFfilename or \"None\"> <parameter name 1> <parameter name 2> ... <parameter value 1> <parameter value 2> ... \n");
      printf("Example: ./main 0 2 None par1 par2 1.0 2.0 \n");
      exit(1);
  }
                               
  //err=readVar(argv[1]);
  flags = strtoul(argv[1],&errstr,10);
  if(*errstr != '\0'){
      printf("invalid input:%s for the number of parameters\n", errstr);
      exit(1);
    }
 
  n_inputvals = strtol(argv[2],&errstr,10);
  if(*errstr != '\0'){
      printf("invalid input:%s for the number of parameters\n", errstr);
      exit(1);
    }
   
 
  if(1 + n_options + 2*n_inputvals != argc){
      printf("invalid input: missmatching the number of parameters (passed: %d, required: 2 * %d)\n", argc-n_options, n_inputvals);
      exit(1);
  }
  /********
   Note: order of parameters:  
      ./main n m fname par1 par2 ... parn val1 val2 ...
      0      1 2 3     4    5    ... n+3  n+4  n+5  ...
  ********/
  err = readVar_str(n_inputvals,&argv[1+n_options],&argv[n_inputvals+1+n_options]);
  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
           
  

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
  if(CDM1) 
  { 
     qNumbers(CDM1, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM1,  spin2,Mcdm1); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  if(CDM2) 
  { 
     qNumbers(CDM2, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM2,spin2,Mcdm2); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }

  
if (flags & MASSES_INFO)
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}


if (flags & CONSTRAINTS)
{ double csLim;
  if(Zinvisible()) printf("Excluded by Z->invizible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q q-\\bar  Cross Section= %.2E pb\n",csLim);
}
 

if (flags & MONOJET)
{ double CL=monoJet();
  printf(" Monojet signal exclusion CL is %.3e\n", CL);   
}  


if ( (flags & HIGGSBOUNDS) + (flags & HIGGSSIGNALS) )
{  int NH0,NHch;  // number of neutral and charged Higgs particles.
    int HB_id[3],HB_result[3];
   double  HB_obsratio[3],HS_observ,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50]; 
   NH0=hbBlocksMO("HB.in",&NHch);
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
if (flags & HIGGSBOUNDS){
   printf("HiggsBounds(%s)\n", HB_version);
   for(int i=0;i<3;i++) printf("  id= %d  result = %d  obsratio=%.2E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
}
if (flags & HIGGSSIGNALS){
   printf("HiggsSignals(%s)\n",HS_version); 
   printf("  Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_observ,HS_chi2, HS_pval);
}
}


if (flags & LILITH)
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char call_lilith[100], Lilith_version[20];
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



if (flags & SMODELS)
{  int result=0;
   double Rvalue=0;
   char analysis[30]={},topology[30]={};
   int LHCrun=LHC8|LHC13;  // LHC8  - 8TeV; LHC13 - 13TeV; LHC8|LHC13 - both  
#include "../include/SMODELS.inc" 
}   
 



if (flags & OMEGA)
{ int fast=1;
  double Beps=1.E-4, cut=0.01;
  double Omega;  
  int i,err; 
 
  if (strcmp(argv[3],"None") != 0){
      printf("load DOF file: %s....\n", argv[3]);
      err = loadHeffGeff(argv[3]);
      if(err<0){
          printf("invalid input: wrong format\n");
          exit(1);
      }
      if(err==0){ 
          printf("invalid input: cannot open %s\n", argv[3]);
          exit(1);
      }
  }
 
  printf("\n==== Calculation of relic density =====\n");
 
  

  if(CDM1 && CDM2) 
  {
  
    Omega= darkOmega2(fast,Beps);
  /*
    displayPlot("vs1120F","T", Tend, Tstart,1,1,"vs1120F", 0,vs1120F,NULL);
    displayPlot("vs2200F","T", Tend, Tstart,1,1,"vs2200F", 0,vs2200F,NULL);
    displayPlot("vs1100F","T", Tend, Tstart,1,1,"vs1100F", 0,vs1100F,NULL);
    displayPlot("vs1210F","T", Tend, Tstart,1,1,"vs1210F", 0,vs1210F,NULL);
    displayPlot("vs1122F","T", Tend, Tstart,1,1,"vs1122F", 0,vs1122F,NULL);
    displayPlot("vs2211F","T", Tend, Tstart,1,1,"vs2211F", 0,vs2211F,NULL);

    displayPlot("vs1110F","T", Tend, Tstart,1,1,"vs1110F", 0,vs1110F,NULL);
    displayPlot("vs2220F","T", Tend, Tstart,1,1,"vs2220F", 0,vs2220F,NULL);
    displayPlot("vs1110F","T", Tend, Tstart,1,1,"vs1110F", 0,vs1112F,NULL);
    displayPlot("vs1222F","T", Tend, Tstart,1,1,"vs1222F", 0,vs1222F,NULL);
    displayPlot("vs1220F","T", Tend, Tstart,1,1,"vs1220F", 0,vs1220F,NULL);
    displayPlot("vs2210F","T", Tend, Tstart,1,1,"vs2210F", 0,vs2210F,NULL);
    displayPlot("vs2221F","T", Tend, Tstart,1,1,"vs2221F", 0,vs2221F,NULL);
    displayPlot("vs1211F","T", Tend, Tstart,1,1,"vs1211F", 0,vs1211F,NULL);
  */
  
    printf("Omega_1h^2=%.2E\n", Omega*(1-fracCDM2));
    printf("Omega_2h^2=%.2E\n", Omega*fracCDM2);
  } else
  {  double Xf;
     Omega=darkOmega(&Xf,fast,Beps,&err);
     printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
     if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
  }
}

 
 

if (flags & FREEZEIN)
{
  double TR=1E10;
  double omegaFi;  
  toFeebleList(CDM1);
  VWdecay=0; VZdecay=0;
  
  omegaFi=darkOmegaFi(TR,&err);
  printf("omega freeze-in=%.3E\n", omegaFi);
  printChannelsFi(0,0,stdout);
}
 



if (flags & INDIRECT_DETECTION)
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



  if(SpA)
  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
if (flags & SHOWPLOTS){
     sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxA);
}
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }

  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
if (flags & SHOWPLOTS){     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxE);
}
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
if (flags & SHOWPLOTS){    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxP);
}
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  


if (flags & RESET_FORMFACTORS)
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


if (flags & CDM_NUCLEON)
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

  if(CDM1)
  {  
    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes for %s \n",CDM1);
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);
  }
  if(CDM2)
  {
    nucleonAmplitudes(CDM2, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes for %s \n",CDM2);
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);
  }     
}

  
if (flags & CDM_NUCLEUS)
{ double dNdE[300];
  double nEvents;

printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
if (flags & SHOWPLOTS){
    displayPlot("Distribution of recoil energy of 73Ge","E[KeV]",1,50,0,1,"dN/dE",0,dNdERecoil,dNdE);
}

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
if (flags & SHOWPLOTS){
    displayPlot("Distribution of recoil energy of 131Xe","E[KeV]",1,50,0,1,"dN/dE",0,dNdERecoil,dNdE);
}

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
if (flags & SHOWPLOTS){
    displayPlot("Distribution of recoil energy of 23Na","E[KeV]",1,50,0,1,"dN/dE",0,dNdERecoil,dNdE);
}

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
if (flags & SHOWPLOTS){
  displayPlot("Distribution of recoil energy of 127I","E[KeV]",1,50,0,1,"dN/dE",0,dNdERecoil,dNdE);
}
  
}


if (flags & NEUTRINO){
if(!CDM1 || !CDM2)
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;
  
 printf("\n===============Neutrino Telescope=======  for  "); 
 if(forSun) printf("Sun\n"); else printf("Earth\n");  

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
if (flags & SHOWPLOTS){
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0, 2,"dnu/dE",0,SpectdNdE,nu,"dnu_bar/dE",0,SpectdNdE,nu_bar);
}
{ 
    printf(" E>%.1E GeV neutrino flux       %.2E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,nu_bar,NULL));  
} 
  
/* Upward events */
  
  muonUpward(nu,nu_bar, mu);
if (flags & SHOWPLOTS){
  displayPlot("Upward muons[1/Year/km^2/GeV]","E",Emin,Mcdm/2, 0,1,"mu",0,SpectdNdE,mu);  
}
    printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));
  
/* Contained events */
  muonContained(nu,nu_bar,1., mu);
if (flags & SHOWPLOTS){ 
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E",Emin,Mcdm,0,1,"",0,SpectdNdE,mu); 
}
  printf(" E>%.1E GeV Contained muon flux %.2E [1/Year/km^3]\n",Emin,spectrInfo(Emin/Mcdm,mu,NULL));
}        
}


if (flags & DECAYS)
{ char*  pname = pdg2name(25);
  txtList L;
  double width; 
  if(pname)
  { 
    width=pWidth(pname,&L);  
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  } 
  
  pname = pdg2name(24);  
  if(pname)
  { 
    width=pWidth(pname,&L);  
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  } 
}            


if (flags & CROSS_SECTIONS)
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
 
 

if (flags & CLEAN){
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
  system("rm -f   smodels.in  smodels.log  smodels.out  summary.*");  
}



  killPlots();
  return 0;
}
