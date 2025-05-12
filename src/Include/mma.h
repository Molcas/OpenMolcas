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
#else
#ifndef ADD_
#define allocmem  allocmem_
#define c_getmem  c_getmem_
#define cptr2woff cptr2woff_
#define woff2cptr woff2cptr_
#define mma_avmem mma_avmem_
#define allomblck allomblck_
#define freemblck freemblck_
#endif
#endif

void *woff2cptr(INT offset);
void *allomblck(char *name,  INT *len);
INT   cptr2woff(void *cptr);
INT   freemblck(void *mblck);
INT   allocmem(INT *size);
INT   c_getmem(char *name, char* Op, char *dtyp, INT *offset, INT *len);
INT   mma_avmem(void);
