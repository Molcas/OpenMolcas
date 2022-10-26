!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE JACORD(HH,EIGVEC,NVEC,NDIM)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION HH(*), EIGVEC(NDIM,NVEC)

      ThrZ=1.0D-14
      DO 100 I=1,NVEC-1
        II=(I*(I+1))/2
        EI=HH(II)
        EMIN=EI
        IMIN=I
        DO 10 J=I+1,NVEC
          JJ=(J*(J+1))/2
          EJ=HH(JJ)
          IF(EJ.GE.EMIN.or.ABS(EJ-EMIN).lt.ThrZ) GOTO 10
          EMIN=EJ
          IMIN=J
  10    CONTINUE
        IF(IMIN.EQ.I) GOTO 100
        HH(II)=EMIN
        JJ=(IMIN*(IMIN+1))/2
        HH(JJ)=EI
        DO 20 K=1,NDIM
          SWAP=EIGVEC(K,I)
          EIGVEC(K,I)=EIGVEC(K,IMIN)
          EIGVEC(K,IMIN)=SWAP
  20    CONTINUE
 100  CONTINUE
      RETURN
      END
