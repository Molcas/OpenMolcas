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
      SUBROUTINE VECSUM(       C,       A,       B,    FACA,    FACB,
     &                      NDIM)
C
C     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),C(*)
*
      IF(FACA.NE.0.0D0.AND.FACB.NE.0.0D0) THEN
        DO 100 I=1,NDIM
          S=FACA*A(I)+FACB*B(I)
          C(I)=S
  100   CONTINUE
*
      ELSE IF(FACA.EQ.0.0D0.AND.FACB.NE.0.0D0) THEN
        DO 200 I=1,NDIM
          S=FACB*B(I)
          C(I)=S
  200   CONTINUE
*
      ELSE IF(FACA.NE.0.0D0.AND.FACB.EQ.0.0D0) THEN
        DO 300 I=1,NDIM
          S=FACA*A(I)
          C(I)=S
  300   CONTINUE
*
      ELSE IF(FACA.EQ.0.0D0.AND.FACB.EQ.0.0D0) THEN
        DO 400 I=1,NDIM
          C(I)=0.0D0
  400   CONTINUE
      END IF
C
      RETURN
      END
