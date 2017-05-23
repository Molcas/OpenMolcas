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
      Subroutine VECP(P1,P2,P3,DNORM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P1(3),P2(3),P3(3)
C
C     Esegue il prodotto vettoriale P3 = P1 x P2
C
      P3(1) = P1(2)*P2(3) - P1(3)*P2(2)
      P3(2) = P1(3)*P2(1) - P1(1)*P2(3)
      P3(3) = P1(1)*P2(2) - P1(2)*P2(1)
      DNORM3 = sqrt(P3(1)*P3(1) + P3(2)*P3(2) + P3(3)*P3(3))
      RETURN
      END
