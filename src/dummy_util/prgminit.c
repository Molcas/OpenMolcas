#ifndef _HAVE_EXTRA_

#include <molcastype.h>

#ifdef _CAPITALS_
#define prgminitc PRGMINITC
#define prgmtranslatec PRGMTRANSLATEC
#define prgmfree PRGMFREE
#define isinmem ISINMEM
#else
#ifndef ADD_
#define prgminitc prgminitc_
#define prgmtranslatec prgmtranslatec_
#define prgmfree prgmfree_
#define isinmem isinmem_
#endif
#endif

void prgminitc(char *, INT *);
void prgmtranslatec (char *, INT *, char *, INT *, INT*);
void prgmfree(void);
INT isinmem(char *);
void setsubdir(char *);

void prgminitc(char *module, INT *lenmodule) {};

void prgmtranslatec (char *in, INT *lin, char *out, INT *lout, INT* ms) {};

void prgmfree(void) {};

#ifndef _GA_
INT isinmem(char *fname) {};
#endif

void setsubdir(char *sub) {};

#endif
