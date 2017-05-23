************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2001, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ABSTR_TO_ORDSTR( IA_OC, IB_OC,  NAEL,  NBEL,IDET_OC,
     &                           IDET_SP,ISIGN)
*
* An alpha string (IA) and a betastring (IB) is given.
* Combine these two strings to give an determinant with
* orbitals in ascending order. For doubly occupied orbitals
* the alphaorbital is given first.
* The output is given as IDET_OC : Orbital occupation (configuration )
*                        IDET_SP : Spin projections

* The phase required to change IA IB into IDET is computes as ISIGN
*
* Jeppe Olsen, November 2001
*
#include "implicit.fh"
*. Input
      INTEGER IA_OC(NAEL),IB_OC(NBEL)
*. Output
      INTEGER IDET_OC(NAEL+NBEL)
      INTEGER IDET_SP(NAEL+NBEL)
*
      NEXT_AL = 1
      NEXT_BE = 1
      NEXT_EL = 0
      ISIGN = 1
*. Loop over next electron in outputstring
      DO NEXT_EL = 1, NAEL+NBEL
       IF(NEXT_AL.LE.NAEL.AND.NEXT_BE.LE.NBEL) THEN
*
         IF(IA_OC(NEXT_AL).LE.IB_OC(NEXT_BE)) THEN
*. Next electron is alpha electron
           IDET_OC(NEXT_EL) = IA_OC(NEXT_AL)
           IDET_SP(NEXT_EL) = +1
           NEXT_AL = NEXT_AL + 1
         ELSE
*. Next electron is beta electron
           IDET_OC(NEXT_EL) = IB_OC(NEXT_BE)
           IDET_SP(NEXT_EL) = -1
           NEXT_BE = NEXT_BE + 1
           ISIGN = ISIGN*(-1)**(NAEL-NEXT_AL+1)
         END IF
       ELSE IF(NEXT_BE.GT.NBEL) THEN
*. Next electron is alpha electron
           IDET_OC(NEXT_EL) = IA_OC(NEXT_AL)
           IDET_SP(NEXT_EL) = +1
           NEXT_AL = NEXT_AL + 1
       ELSE IF(NEXT_AL.GT.NAEL) THEN
*. Next electron is beta electron
           IDET_OC(NEXT_EL) = IB_OC(NEXT_BE)
           IDET_SP(NEXT_EL) = -1
           NEXT_BE = NEXT_BE + 1
           ISIGN = ISIGN*(-1)**(NAEL-NEXT_AL+1)
       END IF
      END DO
*    ^ End of loop over electrons in outputlist
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ABSTR to ORDSTR : '
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' Input alpha and beta strings '
        CALL IWRTMA(IA_OC,1,NAEL,1,NAEL)
        CALL IWRTMA(IB_OC,1,NBEL,1,NBEL)
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(IDET_OC,1,NAEL+NBEL,1,NAEL+NBEL)
        WRITE(6,*) ' Spin projections '
        CALL IWRTMA(IDET_SP,1,NAEL+NBEL,1,NAEL+NBEL)
      END IF
*
      RETURN
      END
