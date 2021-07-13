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
      REAL*8 FUNCTION FNDMNX(VECTOR,NDIM,MINMAX)
!
!     FIND SMALLEST(MINMAX=1) OR LARGEST(MINMAX=2)
!     ABSOLUTE VALUE OF ELEMENTS IN VECTOR
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VECTOR(*)
!
! jwk-cleanup
      RESULT = 0.0D0
      IF(MINMAX.EQ.1) THEN
       RESULT=ABS(VECTOR(1))
       DO 100 I=2,NDIM
        RESULT=MIN(RESULT,ABS(VECTOR(I)))
  100  CONTINUE
      END IF
!
      IF(MINMAX.EQ.2) THEN
       RESULT=ABS(VECTOR(1))
       DO 200 I=2,NDIM
        RESULT=MAX(RESULT,ABS(VECTOR(I)))
  200  CONTINUE
       END IF
!
       FNDMNX=RESULT
      RETURN
      END
