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
      SUBROUTINE MATINV (A,B,N,L,IDIM)
C      IF L=0 RETURNS INVERSE OF A IN A, IF L=1 SOLUTION OF AX=B IN B
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 A(IDIM,IDIM), B(IDIM)
*
*...  internal variables
      PARAMETER (maxdim=44)
      Integer IP(maxdim),IN(maxdim,2)
*
      IF (IDIM.GT.maxdim) THEN
        WRITE (6,*) 'MATINV: Idim',Idim
        WRITE (6,*) 'Abend: Increase maxdim !!'
        Call Abend
      ENDIF
      IR=0
      IC=0
      D=1.D0
      DO  I=1,N
         IP(I)=0
      End Do
      DO 12 I=1,N
      AMAX=0.D0
      DO 3 J=1,N
      IF(IP(J).GT.0) GO TO 3
      IF(IP(J).LT.0) GO TO 4
      DO 2 K=1,N
      IF(IP(K).EQ.1) GO TO 2
      IF(IP(K).GT.1) GO TO 4
      IF(   abs( A(J,K) ) .LE.  AMAX )  GO TO 2
      IR=J
      IC=K
      AMAX =   abs(  A(J,K)  )
    2 CONTINUE
    3 CONTINUE
      IP(IC)=IP(IC)+1
      IF( AMAX.GT.1D-30 ) GO TO 6
    4 WRITE(6,105)
  105 FORMAT(' * '/' * ',16H SINGULAR MATRIX)
      Call Abend
    6 IF(IR.EQ.IC)  GO TO 8
      D=-D
      DO 7 K=1,N
      AMAX=A(IR,K)
      A(IR,K)=A(IC,K)
    7 A(IC,K)=AMAX
      IF(L.EQ.0) GO TO 8
      AMAX=B(IR)
      B(IR)=B(IC)
      B(IC)=AMAX
    8 IN(I,1)=IR
      IN(I,2)=IC
      AMAX=A(IC,IC)
      D=D*AMAX
      A(IC,IC)=1.D0
      DO 9 K=1,N
    9 A(IC,K)=A(IC,K)/AMAX
      IF(L.EQ.0) GO TO 10
      B(IC)=B(IC)/AMAX
   10 DO 12 J=1,N
      IF(J.EQ.IC) GO TO 12
      AMAX=A(J,IC)
      A(J,IC)=0.D0
      DO 11 K=1,N
   11 A(J,K)=A(J,K)-A(IC,K)*AMAX
      IF(L.EQ.0) GO TO 12
      B(J)=B(J)-B(IC)*AMAX
   12 CONTINUE
      IF(L.EQ.1) GO TO 15
      DO 14 I=1,N
      J=N+1-I
      IF(IN(J,1).EQ.IN(J,2)) GO TO 14
      IR=IN(J,1)
      IC=IN(J,2)
      DO 13 K=1,N
      AMAX=A(K,IR)
      A(K,IR)=A(K,IC)
   13 A(K,IC)=AMAX
   14 CONTINUE
   15 CONTINUE
      RETURN
      END
