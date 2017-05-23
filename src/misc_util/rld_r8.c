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
* Copyright (C) 2002, Per-Olof Widmark                                 *
***********************************************************************/
/**************************************************************************/
/*                                                                        */
/* This routine performs run length decoding of floating point numbers    */
/* previously encoded by lre_r8.                                          */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/*          Lund University, Sweden                                       */
/* Written: October 2002                                                  */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <molcastype.h>

#ifdef _CAPITALS_
#define rld_r8 RLD_R8
#else
#ifndef ADD_
#define rld_r8 rld_r8_
#endif
#endif

#ifdef _BIG_ENDIAN_
#define IND 0
#else
#define IND 3
#endif

void rld_r8(double in[], INT *n_in, double out[], INT *n_out, double *thr) {
   unsigned long long int  *ptr_64;
   unsigned short          *ptr_16;
   INT      i,k,m,n,r,cnt;

   if(0) {
      n=0;
      ptr_64=(unsigned long long int *)in;
      ptr_16=(unsigned short *)in;
      for(k=0; k<*n_in; k++) {
         if(ptr_16[4*k+IND]==0) {
            m=ptr_64[k];
            for(i=0; i<m; i++) out[n++]=0.0;
         } else {
            out[n++]=in[k];
         }
      }
      *n_out=n;
   } else {
      n=0;
      cnt=0;
      ptr_64=(unsigned long long int *)in;
      ptr_16=(unsigned short *)in;
      for(k=0; n<*n_out ; k++) {
         if(ptr_16[4*k+IND]==0) {
            m=ptr_64[k];
            r=m;
            for(i=0; (i<m) & (n<*n_out); i++) { out[n++]=0.0; r--; }
            if(r!=0) {
               ptr_64[k]=r;
               cnt--;
            }
         } else {
            out[n++]=in[k];
         }
         cnt++;
      }
      if(0) *n_out=n;
      if(1) *n_in=cnt;
   }

}
