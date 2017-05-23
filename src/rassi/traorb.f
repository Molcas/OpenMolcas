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
      SUBROUTINE TRAORB(NSYM,NOSH,NBASF,NCXA,CXA,NCMO,CMO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NOSH(NSYM),NBASF(NSYM),CXA(NCXA),CMO(NCMO)
C TRANSFORM ORBITAL COEFFICIENTS CMO BY MULTIPLYING WITH
C TRANSFORMATION MATRIX CXA.
#include "WrkSpc.fh"
      CALL GETMEM('CNEW  ','ALLO','REAL',LCNEW,NCMO)
      ISTA1=1
      ISTA2=1
      DO 10 IS=1,NSYM
        NO=NOSH(IS)
        IF(NO.EQ.0) GOTO 10
        NB=NBASF(IS)
        IF(NB.EQ.0) GOTO  5
*        CALL MXMA (CMO(ISTA1),1,NB,CXA(ISTA2),1,NO,
*     *             WORK(LCNEW-1+ISTA1),1,NB,NB,NO,NO)
        CALL DGEMM_('N','N',NB,NO,NO,1.0D0,
     &             CMO(ISTA1),NB,CXA(ISTA2),NO,
     &        0.0D0, WORK(LCNEW-1+ISTA1),NB)
        ISTA1=ISTA1+NO*NB
5       ISTA2=ISTA2+NO**2
10    CONTINUE
      CALL DCOPY_(NCMO,WORK(LCNEW),1,CMO,1)
      CALL GETMEM('      ','FREE','REAL',LCNEW,NCMO)
      RETURN
      END
