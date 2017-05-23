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
/**************************************************************************/
/*                                                                        */
/* This routine figures out how much one ULP is worth.                    */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>

#ifdef _CAPITALS_
#define one_ulp ONE_ULP
#else
#ifndef ADD_
#define one_ulp one_ulp_
#endif
#endif

#ifdef _BIG_ENDIAN_
#define IND 7
#else
#define IND 0
#endif

void one_ulp(double *ulp) {
   union {
      double                 x;
      unsigned char          i8[8];
   } share;
   double ulp1,ulp2,ulp3,ulp4,big_ulp;

   share.x=1.0; share.i8[IND]=1; ulp1=(share.x-1.0)/1.0;
   share.x=2.0; share.i8[IND]=1; ulp2=(share.x-2.0)/2.0;
   share.x=4.0; share.i8[IND]=1; ulp3=(share.x-4.0)/4.0;
   share.x=8.0; share.i8[IND]=1; ulp4=(share.x-8.0)/8.0;
   if(0) printf("ulp1=%6.3e, ulp2=%6.3e, ulp3=%6.3e, ulp4=%6.3e\n",ulp1,ulp2,ulp3,ulp4);

   big_ulp=0.0;
   if(ulp1>big_ulp) big_ulp=ulp1;
   if(ulp2>big_ulp) big_ulp=ulp2;
   if(ulp3>big_ulp) big_ulp=ulp3;
   if(ulp4>big_ulp) big_ulp=ulp4;

   *ulp=big_ulp;

}
