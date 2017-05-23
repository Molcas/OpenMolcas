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
/* This routine counts how long each value is when using the tce          */
/* compression routines.                                                  */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/*          Lund University, Sweden                                       */
/* Written: October 2002                                                  */
/*                                                                        */
/**************************************************************************/
#include <math.h>
#include <molcastype.h>

#ifdef _CAPITALS_
#define tcl_r8 TCL_R8
#define one_ulp ONE_ULP
#else
#ifndef ADD_
#define tcl_r8 tcl_r8_
#define one_ulp one_ulp_
#endif
#endif

#ifdef _BIG_ENDIAN_
#define IND 0
#else
#define IND 3
#endif

void one_ulp(double *);
void tcl_r8(double in[], INT *n_in, INT out[], double *thr, INT *Init_do_setup_l) {
/* unsigned long long int *ptr_64; */
   unsigned short         *ptr_16;
/* unsigned char          *ptr_8; */
/* marked as volatile to prevent over-optimization */
   volatile union {
      double dbl;
      short  i16[4];
   } share;
   static unsigned char tab[65536];
   double ulp;
   static int do_setup=1;
   INT i,k,n,skip,keep;

   if(*Init_do_setup_l == 1) do_setup=1;
   if(do_setup) {
      do_setup=0;
      one_ulp(&ulp);
      for(k=0; k<65536; k++) tab[k]=8; tab[0]=2;
      for(k=0; k<65536; k++) {
/* ifdef _FP_IEEE_ */
         if(k<16)               continue;
         if(k>64879)            continue;
         if(k>32111 && k<32784) continue;
/* endif */
         share.dbl=0.0;
         share.i16[IND]=(short)k;
         skip=0;
         for(i=0; i<6; i++) {
            if(fabs(ulp*share.dbl*255.0)<*thr) skip++;
            else break;
            share.dbl=share.dbl*256.0;
         }
         if(skip>6) skip=6;
         tab[k]=8-skip;
      }
   }

/* ptr_64=(unsigned long long int *)out; */
   n=0;
   for(k=0; k<*n_in; k++) {
/*    ptr_8=(unsigned char *)&in[k]; */
      ptr_16=(unsigned short *)&in[k];
      keep=tab[ptr_16[IND]];
      out[n++]=keep;
   }

}
