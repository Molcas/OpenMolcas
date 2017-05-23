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
      SUBROUTINE ORDER(C,D,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(N,N),D(N)
      DO 30 I=1,N-1
        IMIN=I
        DMIN=D(I)
        DO 10 J=I+1,N
          IF(D(J).GE.DMIN) GOTO 10
          DMIN=D(J)
          IMIN=J
10      CONTINUE
        IF(I.EQ.IMIN) GOTO 30
        D(IMIN)=D(I)
        D(I)=DMIN
        DO 20 K=1,N
          TMP=C(K,I)
          C(K,I)=C(K,IMIN)
          C(K,IMIN)=TMP
20      CONTINUE
30    CONTINUE
      RETURN
      END
