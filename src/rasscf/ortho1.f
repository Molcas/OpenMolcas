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
      SUBROUTINE ORTHO1(S,U,V,N,M)
C
C     Purpose: Orthonormalize N times N vector set U.
C
C     Called from ORTHO.
C
C     Subroutine calls: ORTHO2.
C
C     ********** IBM-3090 Release 88 10 10 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S(*),U(*),V(*)
      THR=0.2D0
      IBASE=1
C
      DO I=1,M

10     CONTINUE
       CALL ORTHO2(S,U(IBASE),V(IBASE),N)
C     Normalize U and calculate V=S*U.
       XNORM=1.0D0
       JBASE=1
       DO J=1,I-1
        OVL=DDOT_(N,U(IBASE),1,V(JBASE),1)
        OVL1=-OVL
        CALL DAXPY_(N,OVL1,U(JBASE),1,U(IBASE),1)
        XNORM=XNORM-OVL**2
        IF (XNORM.LT.THR) GOTO 10
        JBASE=JBASE+N
       END DO
       CALL ORTHO2(S,U(IBASE),V(IBASE),N)
       IBASE=IBASE+N

      END DO

      RETURN
      END
