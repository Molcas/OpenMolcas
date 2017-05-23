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
C
      SUBROUTINE DECOMP(NN,A,UL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NN,NN),UL(NN,NN),SCALES(200)
      COMMON /IPS/IPS(200)
      N=NN
      IDXPIV=0 ! dummy initialize
C**** INITIALIZE IPS, UL AND SCALES
      DO 5 I=1,N
      IPS(I)=I
      ROWNRM=0.0D00
      DO 2 J=1,N
      UL(I,J)=A(I,J)
c      IF (ROWNRM-ABS(UL(I,J))) 1,2,2
      IF (ROWNRM-ABS(UL(I,J)).ge.0) goto 2
c1     CONTINUE
      ROWNRM=abs(UL(I,J))
2     CONTINUE
c      IF (ROWNRM) 3,4,3
      IF (ROWNRM.eq.0) goto 4
c3     CONTINUE
      SCALES(I)=1.0D00/ROWNRM
      GO TO 5
4     CALL SING(1)
      SCALES(I)=0.0D00
5     CONTINUE
C**** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      NM1=N-1
      DO 17 K=1,NM1
      BIG=0.0D00
      DO 11 I=K,N
      IP=IPS(I)
      SIZE=ABS(UL(IP,K))*SCALES(IP)
c      IF (SIZE-BIG) 11,11,10
      IF (SIZE-BIG.le.0) goto 11
c10    CONTINUE
      BIG=SIZE
      IDXPIV=I
11    CONTINUE
c      IF (BIG) 13,12,13
      IF (BIG.ne.0) goto 13
c12    CONTINUE
      CALL SING(2)
      GOTO 17
c13    IF (IDXPIV-K) 14,15,14
13    IF (IDXPIV.eq.K) goto 15
c14    CONTINUE
      J=IPS(K)
      IPS(K)=IPS(IDXPIV)
      IPS(IDXPIV)=J
15    KP=IPS(K)
      PIVOT=UL(KP,K)
      KP1=K+1
      DO 16 I=KP1,N
      IP=IPS(I)
      EM=-UL(IP,K)/PIVOT
      UL(IP,K)=-EM
      DO 16 J=KP1,N
      UL(IP,J)=UL(IP,J)+EM*UL(KP,J)
16    CONTINUE
17    CONTINUE
      KP=IPS(N)
c      IF (UL(KP,N)) 19,18,19
      IF (UL(KP,N).ne.0) goto 19
c18    CONTINUE
      CALL SING(2)
19    RETURN
      END
