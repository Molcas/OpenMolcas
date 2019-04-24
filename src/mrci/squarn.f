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
      SUBROUTINE SQUARN(A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(N,N)
      IN=2
      DO 10 I=2,N
        CALL DCOPY_(I-1,A(IN),1,B(1,I),1)
        IN=IN+I
10    CONTINUE
      DO 20 I=1,N-1
        CALL VNEG(B(I,I+1),N,B(I+1,I),1,N-I)
20    CONTINUE
      CALL DCOPY_(N,[0.0D00],0,B,N+1)
      RETURN
      END
