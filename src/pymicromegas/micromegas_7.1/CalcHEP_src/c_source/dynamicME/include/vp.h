#ifndef __VP__
#define __VP__

#include "../../../include/VandP.h"

extern int      pTabPos(const char * name);
extern double   pMass(const char * name);
extern int      pNum(const char * name);
extern int      qNumbers(const char *pname, int *spin2, int * charge3, int * cdim);
extern char*    pdg2name(int pdg);
extern REAL*    varAddress(const char * name);
extern char*    antiParticle(const char*name); 

#endif
