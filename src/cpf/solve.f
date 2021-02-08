************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
C
      SUBROUTINE SOLVE(NN,UL,B,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION UL(NN,NN),B(*),X(*)
      COMMON /IPS/IPS(200)
      N=NN
      NP1=N+1
C
      IP=IPS(1)
      X(1)=B(IP)
      DO 2 I=2,N
       IP=IPS(I)
       IM1=I-1
       SUM=0.0D00
       DO 1 J=1,IM1
        SUM=SUM+UL(IP,J)*X(J)
1      CONTINUE
       X(I)=B(IP)-SUM
2     CONTINUE
      IP=IPS(N)
      X(N)=X(N)/UL(IP,N)
      DO 4 IBACK=2,N
       I=NP1-IBACK
C****  I GOES (N-1),...,1
       IP=IPS(I)
       IP1=I+1
       SUM=0.0D00
       DO 3 J=IP1,N
        SUM=SUM+UL(IP,J)*X(J)
3      CONTINUE
       X(I)=(X(I)-SUM)/UL(IP,I)
4     CONTINUE
      RETURN
      END
