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
      SUBROUTINE JACSCF (A,B,C,NAA,NQQ,EPSLON)
      IMPLICIT real*8 (A-H,O-Z)
C VERSION 4 AUGUST 1971
C SUBROUTINE TO FIND ALL THE EIGENVALUES AND EIGENVECTORS OF A
C BY THE JACOBI METHOD
C THIS PROGRAM TREATS THE ORTHOGONAL CASE (S=1)
C NAA=DIMENSION OF A,B,C
C A TRIANGULAR, B MATRIX OF VECTORS, C EIGENVALUES
C B CLEARED TO UNIT MATRIX IF NQ=-1, NOT CLEARED IF NQ= NO. ROWS B
C NQ MUST NOT BE LESS THAN NAA
C EPSLON IS THE CONVERGENCE CRITERIA FOR OFF DIAGONAL ELEMENTS
      DIMENSION A(*),B(*),C(*)
*------
* POW: Unnecessary but warning stopping initialization
      term=1.0d30
*------
      NQ=NQQ
      LOOPC=0
      NA=NAA
      NN=(NA*(NA+1))/2
c      IF (NQ) 10,10,50
      IF (NQ.gt.0) goto 50
c10    CONTINUE
      K=1
      NQ=NA
      DO 40 I=1,NA
      DO 41 J=1,NA
c      IF (I-J) 20,30,20
      IF (I.eq.J) goto 30
c20    CONTINUE
      B(K)=0.0D00
      GO TO 42
30    B(K)=1.0D00
42    K=K+1
41    CONTINUE
40    CONTINUE
50    SUM=0.0D00
c      IF (NA-1) 330,310,60
      IF (NA-1.lt.0) goto 330
      IF (NA-1.eq.0) goto 310
c60    CONTINUE
      K=1
      AMAX=0.0D00
      DO 110 I=1,NA
      DO 100 J=1,I
c      IF (I-J) 70,90,70
      IF (I.eq.J) goto 90
c70    IF (ABS(A(K))-AMAX) 90,90,80
      IF (ABS(A(K))-AMAX.le.0) goto 90
c80    CONTINUE
      AMAX=ABS(A(K))
90    TERM=A(K)*A(K)
      SUM=SUM+TERM+TERM
      K=K+1
100   CONTINUE
      SUM=SUM-TERM
110   CONTINUE
      SUM=SQRT(SUM)
      THRESH=SUM/SQRT(dble(NA))
      THRSHG=THRESH*EPSLON
c      IF (THRSHG-AMAX) 120,310,310
      IF (THRSHG-AMAX.ge.0) goto 310
c120   CONTINUE
      THRESH=AMAX/3.
c      IF (THRESH-THRSHG) 130,140,140
      IF (THRESH-THRSHG.ge.0) goto 140
c130   CONTINUE
      THRESH=THRSHG
140   K=2
      N=0
      JD=1
      KDX=0
      DO 250 J=2,NA
      ID=0
      JD=JD+J
      JJ=J-1
      KC=0
      KDX=KDX+NQ
      DO 240 I=1,JJ
      ID=ID+I
      IF (ABS(A(K)).GT.THRESH) GO TO 150
      KC=KC+NQ
      GO TO 241
150   N=N+1
      ALPHA=(A(JD)-A(ID))/(2.*A(K))
      BETA=1.0D00/(1.0D00+ALPHA*ALPHA)
      ROOT=1.0D00+ABS(ALPHA)*SQRT(BETA)
      SSQ=BETA/(2*ROOT)
      CSQ=ROOT/2
      CC=SQRT(CSQ)
      S=-SQRT(SSQ)*SIGN(1.0D00,ALPHA)
      TWOSC=2*CC*S
      TEMPA=CSQ*A(ID)+TWOSC*A(K)+SSQ*A(JD)
      A(JD)=CSQ*A(JD)-TWOSC*A(K)+SSQ*A(ID)
      A(ID)=TEMPA
      A(K)=0.0D00
      KA=JD-J
      KB=ID-I
      KD=KDX
      DO 230 L=1,NQ
      KC=KC+1
      KD=KD+1
      TEMPA=CC*B(KC)+S*B(KD)
      B(KD)=-S*B(KC)+CC*B(KD)
      B(KC)=TEMPA
      IF (L.GT.NA) GO TO 230
c      IF (I-L) 180,160,200
      IF (I-L.lt.0) goto 180
      if (I-L.gt.0) goto 200
c160   CONTINUE
      KB=KB+1
170   KA=KA+1
      GO TO 230
180   KB=KB+L-1
c      IF (J-L) 190,170,210
      IF (J-L.eq.0) goto 170
      IF (J-L.gt.0) goto 210
c190   CONTINUE
      KA=KA+L-1
      GO TO 220
200   KB=KB+1
210   KA=KA+1
220   TEMPA=CC*A(KB)+S*A(KA)
      A(KA)=-S*A(KB)+CC*A(KA)
      A(KB)=TEMPA
230   CONTINUE
241   K=K+1
240   CONTINUE
      K=K+1
250   CONTINUE
      LOOPC=LOOPC+1
c      IF (LOOPC-50) 260,340,340
      IF (LOOPC-50.ge.0) goto 340
c260   IF (N-NN/8) 270,270,140
      IF (N-NN/8.gt.0) goto 140
c270   IF (THRESH-THRSHG) 280,300,280
      IF (THRESH-THRSHG.eq.0) goto 300
c280   CONTINUE
      THRESH=THRESH/3.
c      IF (THRESH-THRSHG) 290,140,140
      IF (THRESH-THRSHG.ge.0) goto 140
c290   CONTINUE
      THRESH=THRSHG
      GO TO 140
c300   IF (N) 140,310,140
300   IF (N.ne.0) goto 140
310   LL=0
      DO 320 L=1,NA
      LL=LL+L
      C(L)=A(LL)
320   CONTINUE
340   CONTINUE
330   RETURN
      END
