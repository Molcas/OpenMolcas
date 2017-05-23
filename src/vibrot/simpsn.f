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
      SUBROUTINE SIMPSN(X,H,NDIM,F)
C
C     THIS SUBROUTINE EVALUATES THE INTEGRAL INT(XDX) USING SIMPSONS
C     RULE. THE NUMBER OF STEPS IS NDIM (ASSUMED ODD) WITH THE GRID
C     SIZE H. THE RESULTING VALUE OF THE INTEGRAL IS RETURNED IN F.
C
C     ********** RELEASE 81 04 09 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*)
      DATA D2,D6/2.D 00,6.D 00/
C
      F=X(1)+X(NDIM)
      NDIM1=NDIM-1
      DO 10 I=2,NDIM1
      F=F+D2*X(I)
      IF(MOD(I,2).EQ.0) F=F+D2*X(I)
10    CONTINUE
      F=F*H/D6
      RETURN
      END
