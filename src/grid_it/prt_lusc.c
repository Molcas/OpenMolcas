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

/***********************************************************************
* Adapted from SAGIT to work with OpenMolcas (October 2020)            *
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <molcastype.h>

#ifdef _CAPITALS_
#define prt_lusc PRT_LUSC
#define lusopen LUSOPEN
#define dump_lusc DUMP_LUSC
#define prgmtranslatec PRGMTRANSLATEC
#else
#ifndef ADD_
#define prt_lusc prt_lusc_
#define lusopen lusopen_
#define dump_lusc dump_lusc_
#define prgmtranslatec prgmtranslatec_
#endif
#endif

#define MYMAXPATH 1024

void prgmtranslatec (char *, INT *, char *, INT *, INT*);

INT prt_lusc_(INT *lid, char *line, INT *len, INT *isBin)
{
  INT marker[1];
  FILE *fp;
  int i;
  fp = (FILE*) *lid;

// printf("here %ld %ld \n", *len, *isBin );

  if(*isBin==1){
  marker[0]=*len;
  fwrite(marker,sizeof(INT),1, fp);
  }


  for (i = 0; i < *len; i++){ fputc(line[i], fp);}
  if(*isBin==1){
  marker[0]=*len;
  fwrite(marker,sizeof(INT),1, fp);
  }
  else
  {
  fputc(10, fp);
  }
  return 0;
}

INT lusopen(INT *lid, char *fname, INT *fname_len)
{
  FILE *fp;
  char my[MYMAXPATH];
  int slash = '/';
  char *ptr;
  INT tmp0, tmp1;
  INT ms = 1;
  fname[*fname_len] = 0;

  tmp0 = strlen(fname);
  ptr = strchr(fname, slash);
  if (ptr != NULL) {
    strncpy(my,fname,MYMAXPATH-1);
    my[*fname_len] = 0;
  } else {
    prgmtranslatec(fname, &tmp0, my, &tmp1, &ms);
  }

  fp = fopen(my, "wb");
  *lid = (INT) fp;
  if (fp != NULL) {
    return 0;
  } else {
    return 1;
  }
}
void dump_lusc_(INT *lid, double *Buff, INT *nBuff)
{
 FILE *fp;
  fp = (FILE*) *lid;
  fwrite(Buff,sizeof(double),*nBuff,fp);

 return;
}
