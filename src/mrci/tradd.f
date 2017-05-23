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
      SUBROUTINE TRADD(A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),B(*)
      IN=0
      DO 10 I=1,N
      DO 20 J=1,I
      IN=IN+1
      B(IN)=B(IN)+A(I,J)-A(J,I)
20    CONTINUE
10    CONTINUE
      RETURN
      END
