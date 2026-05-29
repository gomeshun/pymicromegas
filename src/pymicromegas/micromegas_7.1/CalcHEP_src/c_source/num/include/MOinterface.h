#ifndef __MOinterface_
#define __MOinterface_

#include<stdlib.h>
#include<string.h> 
#include<math.h>

//#include"nType.h"
#include"../../../include/num_out.h"

#include"../../ntools/include/vegas.h"
#include"../../dynamicME/include/dynamic_cs.h"

extern int clearHists(void);

extern void  delCuts(void);
extern void delHistograms(void);

extern int  addCut(int key,char*param, double min, double max);
extern int  addHistogram(char*param, double min, double max);
extern int  saveHistograms(char *fname);
//extern void  cutInit(void);
extern vegasGrid * initMCsession( numout*cc, int nSubProc, double P1,double P2, char*STF1, char*STF2);

#endif