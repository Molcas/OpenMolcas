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
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
      Integer(kind=iwp), intent(in):: NSYM, NCXA, NCMO
      Integer(kind=iwp), intent(in):: NOSH(NSYM),NBASF(NSYM)
      Real(kind=wp), intent(in):: CXA(NCXA)
      Real(kind=wp), intent(inout):: CMO(NCMO)

      Real(kind=wp), Allocatable:: CNew(:)
      Integer(kind=iwp) ISTA1,ISTA2,IS,NO,NB

C TRANSFORM ORBITAL COEFFICIENTS CMO BY MULTIPLYING WITH
C TRANSFORMATION MATRIX CXA.
      CALL mma_allocate(CNEW,NCMO,Label='CNEW')
      ISTA1=1
      ISTA2=1
      DO IS=1,NSYM
        NO=NOSH(IS)
        IF(NO.EQ.0) CYCLE
        NB=NBASF(IS)
        IF(NB/=0) THEN
        CALL DGEMM_('N','N',NB,NO,NO,
     &              One,CMO(ISTA1),NB,
     &                    CXA(ISTA2),NO,
     &              Zero,CNEW(ISTA1),NB)
        ISTA1=ISTA1+NO*NB
        END IF
        ISTA2=ISTA2+NO**2
      END DO
      CALL DCOPY_(NCMO,CNEW,1,CMO,1)
      CALL mma_deallocate(CNEW)
      END SUBROUTINE TRAORB
