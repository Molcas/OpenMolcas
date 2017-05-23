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
/* This routine performs run length encoding of floating point numbers    */
/* where all number less than a given threshold are treated as zeros.     */
/* These zeros are run length encoded.                                    */
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
#define rle_r8 RLE_R8
#define one_ulp ONE_ULP
#else
#ifndef ADD_
#define rle_r8 rle_r8_
#define one_ulp one_ulp_
#endif
#endif

#ifdef _BIG_ENDIAN_
#define IND 0
#else
#define IND 3
#endif

void one_ulp(double *);
void rle_r8(double in[], INT *n_in, double out[], INT *n_out, double *thr) {
   unsigned long long int *ptr_64;
   unsigned short         *ptr_16;
   unsigned long long int mask[9]={0xffffffffffffffffLL,
                                   0xffffffffffffff00LL,
                                   0xffffffffffff0000LL,
                                   0xffffffffff000000LL,
                                   0xffffffff00000000LL,
                                   0xffffff0000000000LL,
                                   0xffff000000000000LL,
                                   0xffff000000000000LL,
                                   0xffff000000000000LL};
/* marked as volatile to prevent over-optimization */
   volatile union {
      double dbl;
      short  i16[4];
   } share;
   static unsigned char tab[65536];
   double ulp;
   double screen[9];
   unsigned long long int count;
   static int do_setup=1;
   INT i,k,n,skip;

   if(0) {
      n=0;
      count=0;
      for(k=0; k<*n_in; k++) {
         if(fabs(in[k])>*thr) {
            if(count>0) { ptr_64[n++]=count; count=0; }
            out[n++]=in[k];
            continue;
         }
         count++;
      }
   } else if(0) {
      if(do_setup) {
         do_setup=0;
         one_ulp(&ulp);
         screen[0]=*thr/ulp/255.0;
         for(k=1; k<9; k++) screen[k]=screen[k-1]/256.0;
      }
      ptr_64=(unsigned long long int *)out;
      n=0;
      count=0;
      for(k=0; k<*n_in; k++) {
         if(fabs(in[k])>*thr) {
            if(count>0) { ptr_64[n++]=count; count=0; }
            if(fabs(in[k]) > screen[0])      skip=0;
            else if(fabs(in[k]) > screen[1]) skip=1;
            else if(fabs(in[k]) > screen[2]) skip=2;
            else if(fabs(in[k]) > screen[3]) skip=3;
            else if(fabs(in[k]) > screen[4]) skip=4;
            else if(fabs(in[k]) > screen[5]) skip=5;
            else                             skip=6;
            out[n]=in[k];
            ptr_64[n]=(mask[skip]&ptr_64[n]);
            n++;
            continue;
         }
         count++;
      }
   } else {
      if(do_setup) {
         do_setup=0;
         one_ulp(&ulp);
         for(k=0; k<65536; k++) tab[k]=0;
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
            tab[k]=skip;
         }
      }
      ptr_64=(unsigned long long int *)out;
      n=0;
      count=0;
      for(k=0; k<*n_in; k++) {
         if(fabs(in[k])>*thr) {
            if(count>0) { ptr_64[n++]=count; count=0; }
            ptr_16=(unsigned short *)&in[k];
            skip=tab[ptr_16[IND]];
            out[n]=in[k];
/*
            ptr_64[n]=(mask[skip]&ptr_64[n]);
*/
            n++;
            continue;
         }
         count++;
      }
   }
   if(count>0) { ptr_64[n++]=count; count=0; }
   *n_out=n;

}
