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
      SUBROUTINE FILLMA(N,S,OVE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S(2),OVE(N,N)
      IJ=0
      DO 1 I=1,N
      DO 2 J=1,I
      IJ=IJ+1
      OVE(I,J)=S(IJ)
2     OVE(J,I)=S(IJ)
1     CONTINUE
      RETURN
      END
