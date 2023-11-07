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
! NOTE: THE FOLLOWING ROUTINE MAY NOT BE VECTORIZED.
! THERE IS A COMPILER BUG ON FORTRAN VERSION 2.2.0 (JUNE 1987).
! WRONG RESULTS PRODUCED EVEN ON VECTORIZED LEVEL 1.
SUBROUTINE MKCXAL(NDIMEN,TRAL,CXAL)
use Constants, only: Zero, One
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION TRAL(NDIMEN,NDIMEN),CXAL(NDIMEN,NDIMEN)
DO I=1,NDIMEN
  DO J=I,NDIMEN
    CXAL(I,J)=Zero
  END DO
  CXAL(I,I)=One
END DO
DO K=1,NDIMEN
  DO I=1,K-1
    SUMMA=Zero
    DO J=1,K-1
      SUMMA=SUMMA+CXAL(I,J)*TRAL(J,K)
    END DO
    CXAL(I,K)=-(SUMMA/TRAL(K,K))
  END DO
  DO I=K,NDIMEN
    SUMMA=TRAL(I,K)
    IF (I.EQ.K) SUMMA=-One
    DO J=1,K-1
      SUMMA=SUMMA+CXAL(I,J)*TRAL(J,K)
    END DO
    CXAL(I,K)=-(SUMMA/TRAL(K,K))
  END DO
END DO

END SUBROUTINE MKCXAL
