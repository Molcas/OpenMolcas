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
      !SUBROUTINE SQPRT(A,N)
      !IMPLICIT REAL*8 (A-H,O-Z)
      !DIMENSION A(N,N)
      !CHARACTER*60 FMT
      !BIG=0.0D0
      !DO 10 I=1,N
      !DO 11 J=1,N
      !   BIG=MAX(BIG,ABS(A(I,J)))
!11   !CONTINUE
!10   !CONTINUE
      !IF(0.1D0.LT.BIG .AND. BIG.LT.10000.0D0) THEN
      !   FMT='(8(1X,F12.6))'
      !ELSE
      !   FMT='(8(1X,E12.6))'
      !END IF
      !DO 20 I=1,N
      !   WRITE(6,FMT) (A(I,J),J=1,N)
!20   !CONTINUE
      !RETURN
      SUBROUTINE SQPRT(V,N)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION V(N,N)
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C     -V- IS -N- ROWS BY -M- COLUMNS, WITH LEADING DIMENSION -NDIM-
C
      MAX = 5
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (6,9008)
      WRITE (6,9028) (I,I = IMIN,IMAX)
      WRITE (6,9008)
      DO 120 J = 1,N
  120 WRITE (6,9048) J,(V(J,I),I = IMIN,IMAX)
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(6X,10(4X,I4,4X))
 9048 FORMAT(I5,1X,10F12.7)
      END
