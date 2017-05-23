/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
#ifdef _CAPITALS_
#define allocmem  ALLOCMEM
#define c_getmem  C_GETMEM
#define cptr2woff CPTR2WOFF
#define woff2cptr WOFF2CPTR
#define mma_avmem MMA_AVMEM
#define allomblck ALLOMBLCK
#define freemblck FREEMBLCK
#define lengmblck LENGMBLCK
#define pinnmblck PINNMBLCK
#ifdef _TRACK_
#define trckmblck TRCKMBLCK
#endif
#else
#ifndef ADD_
#define allocmem  allocmem_
#define c_getmem  c_getmem_
#define cptr2woff cptr2woff_
#define woff2cptr woff2cptr_
#define mma_avmem mma_avmem_
#define allomblck allomblck_
#define freemblck freemblck_
#define lengmblck lengmblck_
#define pinnmblck pinnmblck_
#ifdef _TRACK_
#define trckmblck trckmblck_
#endif
#endif
#endif

char *woff2cptr(char etyp[], INT offset);
char *allomblck(char *name,  INT *len);
char *pinnmblck(char *name,  INT *len);
INT   cptr2woff(char etyp[], void *cptr);
INT   freemblck(char *mblck);
INT   lengmblck(char *mblck);
INT   trckmblck(char *mblck);
INT   allocmem(double ref[],char cref[],INT *intof,INT *dblof,INT *sglof, INT *chrof,INT *size);
INT   c_getmem(char *name, char* Op, char *dtyp, INT *offset, INT *len);
INT   mma_avmem(void);
