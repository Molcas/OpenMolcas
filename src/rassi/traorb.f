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
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer NSYM, NCXA, NCMO
      Integer NOSH(NSYM),NBASF(NSYM)
      Real*8 CXA(NCXA),CMO(NCMO)

      Real*8, Allocatable:: CNew(:)

C TRANSFORM ORBITAL COEFFICIENTS CMO BY MULTIPLYING WITH
C TRANSFORMATION MATRIX CXA.
      CALL mma_allocate(CNEW,NCMO,Label='CNEW')
      ISTA1=1
      ISTA2=1
      DO 10 IS=1,NSYM
        NO=NOSH(IS)
        IF(NO.EQ.0) GOTO 10
        NB=NBASF(IS)
        IF(NB.EQ.0) GOTO  5
        CALL DGEMM_('N','N',NB,NO,NO,
     &              1.0D0,CMO(ISTA1),NB,
     &                    CXA(ISTA2),NO,
     &              0.0D0,CNEW(ISTA1),NB)
        ISTA1=ISTA1+NO*NB
5       ISTA2=ISTA2+NO**2
10    CONTINUE
      CALL DCOPY_(NCMO,CNEW,1,CMO,1)
      CALL mma_deallocate(CNEW)
      END SUBROUTINE TRAORB
