/*====== Runtime flags used by pymicromegas Project ======*/

#define MASSES_INFO (1UL << 0)
#define CONSTRAINTS (1UL << 1)
#define MONOJET (1UL << 2)
#define HIGGSBOUNDS (1UL << 3)
#define HIGGSSIGNALS (1UL << 4)
#define LILITH (1UL << 5)
#define SMODELS (1UL << 6)
#define OMEGA (1UL << 7)
#define FREEZEIN (1UL << 8)
#define INDIRECT_DETECTION (1UL << 9)
#define RESET_FORMFACTORS (1UL << 10)
#define CDM_NUCLEON (1UL << 11)
#define CDM_NUCLEUS (1UL << 12)
#define NEUTRINO (1UL << 13)
#define DECAYS (1UL << 14)
#define CROSS_SECTIONS (1UL << 15)
#define SHOWPLOTS (1UL << 16)
#define CLEAN (1UL << 17)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#include "lib/pmodel.h"

int readVar_str(int nVar, char **varnames, char **strvals)
{
  int err = 0;
  char *errstr;
  double val;

  for(int i = 0; i < nVar; i++)
  {
    val = strtod(strvals[i], &errstr);
    if(*errstr != '\0')
    {
      printf("invalid input:%s in %dth argument\n", errstr, i + 1);
      exit(1);
    }
    err = assignVal(varnames[i], val);
    if(err == 1) return i + 1;
  }
  return err;
}

char *pymicromegas_cdm_name(int sector)
{
  if(CDM && sector >= 1 && sector <= Ncdm && CDM[sector]) return CDM[sector];
  if(sector == 1) return CDM1;
  if(sector == 2) return CDM2;
  return NULL;
}

double pymicromegas_cdm_mass(int sector)
{
  if(McdmN && sector >= 1 && sector <= Ncdm) return McdmN[sector];
  if(sector == 1) return Mcdm;
  return 0.0;
}

double pymicromegas_cdm_fraction(int sector)
{
  if(fracCDM && sector >= 1 && sector <= Ncdm) return fracCDM[sector];
  return sector == 1 ? 1.0 : 0.0;
}

void print_cdm_candidate(int sector)
{
  char *name = pymicromegas_cdm_name(sector);
  int spin2, charge3, cdim;

  if(!name) return;
  qNumbers(name, &spin2, &charge3, &cdim);
  printf("\nDark matter candidate %d is '%s' with spin=%d/2 mass=%.2E\n", sector, name, spin2, pymicromegas_cdm_mass(sector));
  if(charge3) printf("Dark Matter has electric charge %d/3\n", charge3);
  if(cdim != 1) printf("Dark Matter is a color particle\n");
}

void load_dof_file(const char *dof_fname)
{
  int err;
  if(strcmp(dof_fname, "None") == 0) return;

  printf("load DOF file: %s....\n", dof_fname);
  err = loadHeffGeff((char *)dof_fname);
  if(err < 0)
  {
    printf("invalid input: wrong format\n");
    exit(1);
  }
  if(err == 0)
  {
    printf("invalid input: cannot open %s\n", dof_fname);
    exit(1);
  }
}

void run_omega(const char *dof_fname)
{
  int err = 0;
  int fast = 1;
  double Beps = 1.E-4;
  double cut = 0.01;
  double Omega;

  load_dof_file(dof_fname);
  printf("\n==== Calculation of relic density =====\n");

  if(Ncdm <= 1)
  {
    double Xf;
    Omega = darkOmega(&Xf, fast, Beps, &err);
    printf("Xf=%.2e Omega=%.2e\n", Xf, Omega);
    if(Omega > 0) printChannels(Xf, cut, Beps, 1, stdout);
  }
  else if(Ncdm == 2)
  {
    Omega = darkOmega2((double)fast, Beps, &err);
    printf("Omega=%.2e err=%d\n", Omega, err);
    for(int sector = 1; sector <= Ncdm; sector++)
    {
      printf("Omega_%dh^2=%.2E\n", sector, Omega * pymicromegas_cdm_fraction(sector));
    }
  }
  else
  {
    Omega = darkOmegaN(fast, Beps, &err);
    printf("Omega=%.2e err=%d\n", Omega, err);
    for(int sector = 1; sector <= Ncdm; sector++)
    {
      printf("Omega_%dh^2=%.2E\n", sector, Omega * pymicromegas_cdm_fraction(sector));
    }
  }
}

void run_freezein(void)
{
  int err = 0;
  double TR = 1E10;
  char *cdm1 = pymicromegas_cdm_name(1);

  if(!cdm1)
  {
    printf("No feeble particle candidate is available\n");
    return;
  }
  toFeebleList(cdm1);
  VWdecay = 0;
  VZdecay = 0;
  printf("omega freeze-in=%.3E\n", darkOmegaFi(TR, cdm1, &err));
  if(err) printf("darkOmegaFi error code = %d\n", err);
  printChannelsFi(0, 0, stdout);
}

void run_indirect_detection(unsigned long flags)
{
  int err;
  double Emin = 1;
  double sigmaV;
  double SpA[NZ], SpE[NZ], SpP[NZ];
  double FluxA[NZ], FluxE[NZ], FluxP[NZ];
  double *SpNe = NULL, *SpNm = NULL, *SpNl = NULL;
  double mass = pymicromegas_cdm_mass(1);

  if(mass > 0) Mcdm = mass;
  printf("\n==== Indirect detection =======\n");
  sigmaV = calcSpectrum(1 + 2 + 4, SpA, SpE, SpP, SpNe, SpNm, SpNl, &err);
  if(err) printf("calcSpectrum error code = %d\n", err);
  if(PlanckCMB(sigmaV, SpA, SpE) > 1) printf("Excluded at 95%% confidence level (1506.03811)\n");

  gammaFluxTab(0.1, 0.05, sigmaV, SpA, FluxA);
  printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n", SpectdNdE(Mcdm / 2, FluxA), Mcdm / 2);
  posiFluxTab(Emin, sigmaV, SpE, FluxE);
  printf("Positron flux = %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV]\n", SpectdNdE(Mcdm / 2, FluxE), Mcdm / 2);
  pbarFluxTab(Emin, sigmaV, SpP, FluxP);
  printf("Antiproton flux = %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV]\n", SpectdNdE(Mcdm / 2, FluxP), Mcdm / 2);
}

void run_cdm_nucleon(void)
{
  char *cdm1 = pymicromegas_cdm_name(1);
  double pA0[2], pA5[2], nA0[2], nA5[2];
  double Nmass = 0.939;
  double mass = pymicromegas_cdm_mass(1);
  double SCcoeff;
  int err;

  if(!cdm1) return;
  printf("\n==== Calculation of CDM-nucleons amplitudes =====\n");
  err = nucleonAmplitudes(cdm1, pA0, pA5, nA0, nA5);
  printf("%s-nucleon micrOMEGAs amplitudes (err=%d)\n", cdm1, err);
  printf("proton:  SI  %.3E  SD  %.3E\n", pA0[0], pA5[0]);
  printf("neutron: SI  %.3E  SD  %.3E\n", nA0[0], nA5[0]);
  if(mass <= 0) return;
  SCcoeff = 4 / M_PI * 3.8937966E8 * pow(Nmass * mass / (Nmass + mass), 2.);
  printf("\n%s-nucleon cross sections[pb] ====\n", cdm1);
  printf(" proton  SI %.3E  SD %.3E\n", SCcoeff * pA0[0] * pA0[0], 3 * SCcoeff * pA5[0] * pA5[0]);
  printf(" neutron SI %.3E  SD %.3E\n", SCcoeff * nA0[0] * nA0[0], 3 * SCcoeff * nA5[0] * nA5[0]);
}

void run_neutrino(void)
{
  double nu[NZ], nu_bar[NZ], mu[NZ];
  double Emin = 1;
  int err;

  if(Ncdm > 1)
  {
    printf("Neutrino telescope output is skipped for multi-component DM.\n");
    return;
  }
  printf("\n===============Neutrino Telescope======= for Sun\n");
  err = neutrinoFlux(Maxwell, 1, nu, nu_bar);
  if(err) printf("neutrinoFlux error code = %d\n", err);
  printf(" E>%.1E GeV neutrino flux       %.2E [1/Year/km^2]\n", Emin, spectrInfo(Emin, nu, NULL));
  printf(" E>%.1E GeV anti-neutrino flux  %.2E [1/Year/km^2]\n", Emin, spectrInfo(Emin, nu_bar, NULL));
  muonUpward(nu, nu_bar, mu);
  printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n", Emin, spectrInfo(Emin, mu, NULL));
}

void run_decay_summary(void)
{
  char *pname = pdg2name(25);
  txtList L;
  double width;

  if(!pname) return;
  width = pWidth(pname, &L);
  printf("\n%s : total width=%E and Branchings:\n", pname, width);
  printTxtList(L, stdout);
}

void run_cross_section_summary(void)
{
  numout *cc;
  int ntot, nin, nout;

  printf("\n====== Example cross sections ======\n");
  cc = newProcess((char *)"e%,E%->2*x");
  if(!cc)
  {
    printf("No e%%,E%%->2*x process is available for this model.\n");
    return;
  }
  procInfo1(cc, &ntot, &nin, &nout);
  printf("Available subprocesses: %d, nin=%d, nout=%d\n", ntot, nin, nout);
}

int main(int argc, char **argv)
{
  int err = 0;
  char *errstr;
  char cdmName[64];
  unsigned long flags;
  int n_inputvals;
  const int n_options = 3;

  ForceUG = 0;
  VZdecay = 0;
  VWdecay = 0;

  if(argc < 1 + n_options)
  {
    printf("usage: ./pymicromegas_main <integer for flags> <the number of parameters> <DOFfilename or None> <parameter names...> <parameter values...>\n");
    printf("Example: ./pymicromegas_main 128 2 None par1 par2 1.0 2.0\n");
    exit(1);
  }

  flags = strtoul(argv[1], &errstr, 10);
  if(*errstr != '\0')
  {
    printf("invalid input:%s for flags\n", errstr);
    exit(1);
  }

  n_inputvals = strtol(argv[2], &errstr, 10);
  if(*errstr != '\0')
  {
    printf("invalid input:%s for the number of parameters\n", errstr);
    exit(1);
  }

  if(1 + n_options + 2 * n_inputvals != argc)
  {
    printf("invalid input: mismatching the number of parameters (passed: %d, required: 2 * %d)\n", argc - 1 - n_options, n_inputvals);
    exit(1);
  }

  err = readVar_str(n_inputvals, &argv[1 + n_options], &argv[n_inputvals + 1 + n_options]);
  if(err > 0)
  {
    printf("Wrong parameter value for %s\n", argv[n_options + err]);
    exit(1);
  }

  err = sortOddParticles(cdmName);
  if(err)
  {
    printf("Can't calculate %s\n", cdmName);
    return 1;
  }

  for(int sector = 1; sector <= Ncdm; sector++) print_cdm_candidate(sector);

  if(flags & MASSES_INFO)
  {
    printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
    printHiggs(stdout);
    printMasses(stdout, 1);
  }

  if(flags & CONSTRAINTS)
  {
    double csLim;
    if(Zinvisible()) printf("Excluded by Z->invisible\n");
    if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q q-bar Cross Section= %.2E pb\n", csLim);
  }

  if(flags & MONOJET)
  {
    printf("Monojet signal exclusion CL is %.3e\n", monoJet());
  }

  if(flags & (HIGGSBOUNDS | HIGGSSIGNALS | LILITH | SMODELS))
  {
    printf("External Higgs/LHC constraint helpers are not run by pymicromegas_main. Use the native model main.c for those workflows.\n");
  }

  if(flags & OMEGA) run_omega(argv[3]);
  if(flags & FREEZEIN) run_freezein();
  if(flags & INDIRECT_DETECTION) run_indirect_detection(flags);

  if(flags & RESET_FORMFACTORS)
  {
    printf("\n======== RESET_FORMFACTORS ======\n");
    calcScalarFF(0.553, 18.9, 70., 35.);
    printf("protonFF d %.2E, u %.2E, s %.2E\n", ScalarFFPd, ScalarFFPu, ScalarFFPs);
    printf("neutronFF d %.2E, u %.2E, s %.2E\n", ScalarFFNd, ScalarFFNu, ScalarFFNs);
  }

  if(flags & CDM_NUCLEON) run_cdm_nucleon();

  if(flags & CDM_NUCLEUS)
  {
    char *expName;
    double pval;
    printf("\n===== Direct detection exclusion ======\n");
    pval = DD_pval(AllDDexp, Maxwell, &expName);
    if(pval < 0.1) printf("Excluded by %s %.1f%%\n", expName, 100 * (1 - pval));
    else printf("Not excluded by DD experiments at 90%% level\n");
  }

  if(flags & NEUTRINO) run_neutrino();
  if(flags & DECAYS) run_decay_summary();
  if(flags & CROSS_SECTIONS) run_cross_section_summary();
  if(flags & CLEAN) cleanDecayTable();

  return 0;
}