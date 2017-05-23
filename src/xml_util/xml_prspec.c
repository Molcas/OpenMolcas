/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Per-Olof Widmark                                       *
***********************************************************************/
/**************************************************************************/
/*                                                                        */
/* This routine prints an xml tag specifier.                              */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/*          Lund University, Sweden                                       */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include "molcastype.h"
#include "xmlapi.h"
#ifdef _CAPITALS_
#define xml_prspec XML_PRSPEC
#else
#ifndef ADD_
#define xml_prcpec xml_prspec_
#endif
#endif
void xml_prspec(FILE *f, char *key, char *value, INT n) {
   char  line[256];
   int   k,m;

   if(n<1) return;
   for(k=0; k<n; k++) line[k]=value[k];
   m=0;
   for(k=0; k<n; k++) if(line[k]!=' ') m=k;
   if(m==0) return;
   line[m+1]=0;
   fprintf(f," %s=\"%s\"",key,line);
}
