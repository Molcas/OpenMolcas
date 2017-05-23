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
      SUBROUTINE VMOV(A,LA,B,LB,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*)
      DO 10 I=0,N-1
         B(1+I*LB)=A(1+I*LA)
10    CONTINUE
      RETURN
      END
