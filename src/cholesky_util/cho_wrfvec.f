************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CHO_WRFVEC(VEC,ISYMA,ISYMB,IVEC1,NUMV)
C
C     Purpose: write full storage vectors to disk.
C
      IMPLICIT NONE
      REAL*8 VEC(*)
      INTEGER ISYMA, ISYMB, IVEC1, NUMV
#include "choreo.fh"

      INTEGER     IOPT, IADR, NTOT

      IOPT = 1
      IADR = NABPK(ISYMA,ISYMB)*(IVEC1 - 1) + 1
      NTOT = NABPK(ISYMA,ISYMB)*NUMV
      CALL DDAFILE(LUFV(ISYMA,ISYMB),IOPT,VEC,NTOT,IADR)

      END
