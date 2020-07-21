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
      SUBROUTINE JACOB_REL(A,B,EIG,N,RANGE,IC)
      IMPLICIT REAL*8(A-H,O-Z)
C     IC=1 EIGNEVALUES AND VECTORS REARRANGED
C       =0 LEFT AS THEY ARE
      DIMENSION A(N,N),B(N,N),EIG(N)
      ENUI=0.D0
      U1=DBLE(N)
      DO 5 I=1,N
      B(I,I)=1.D0
      EIG(I)=A(I,I)
      DO 6 J=1,I
      IF(I.EQ.J)GO TO 6
      B(I,J)=0.D0
      B(J,I)=0.D0
      ENUI=ENUI+A(I,J)*A(I,J)
    6 CONTINUE
    5 CONTINUE
      IF (ENUI.LE.0) THEN
        GO TO 200
      ELSE
        GO TO 10
      END IF
   10 ENUI=sqrt(2.D0*ENUI)
      ENUF=ENUI*RANGE/U1
      IND=0
      THR=ENUI
   15 THR=THR/U1
   20 L=1
   25 M=L+1
   30 IF (abs(A(M,L))-THR.LT.0) THEN
        GO TO 90
      ELSE
        GO TO 35
      END IF
   35 IND=1
      X=.5D0*(EIG(L)-EIG(M))
      Y=-A(M,L)/sqrt(A(M,L)*A(M,L)+X*X)
      IF (X.LT.0) THEN
        GO TO 40
      ELSE
        GO TO 45
      END IF
   40 Y=-Y
   45 IF(Y.GT.1.D0)Y=1.D0
      IF(Y.LT.-1.D0)Y=-1.D0
      XY=1.D0-Y*Y
      SINT=Y/sqrt(2.D0*(1.D0+sqrt(XY)))
      SINT2=SINT*SINT
      COST2=1.D0-SINT2
      COST=sqrt(COST2)
      SINCS=SINT*COST
      DO 85 I=1,N
      IF (I-M.LT.0) THEN
        GO TO 50
      ELSE IF (I-M.EQ.0) THEN
        GO TO 80
      ELSE
        GO TO 55
      END IF
   50 IM=M
      MM=I
      GO TO 60
   55 IM=I
      MM=M
   60 IF (I-L.LT.0) THEN
        GO TO 65
      ELSE IF (I-L.EQ.0) THEN
        GO TO 80
      ELSE
        GO TO 70
      END IF
   65 IL=L
      LL=I
      GO TO 75
   70 IL=I
      LL=L
   75 X=A(IL,LL)*COST-A(IM,MM)*SINT
      A(IM,MM)=A(IL,LL)*SINT+A(IM,MM)*COST
      A(IL,LL)=X
   80 X=B(I,L)*COST-B(I,M)*SINT
      B(I,M)=B(I,L)*SINT+B(I,M)*COST
      B(I,L)=X
   85 CONTINUE
      X=2.D0*A(M,L)*SINCS
      Y=EIG(L)*COST2+EIG(M)*SINT2-X
      X=EIG(L)*SINT2+EIG(M)*COST2+X
      A(M,L)=(EIG(L)-EIG(M))*SINCS+A(M,L)*(COST2-SINT2)
      EIG(L)=Y
      EIG(M)=X
   90 IF (M-N.EQ.0) THEN
        GO TO 100
      ELSE
        GO TO 95
      END IF
   95 M=M+1
      GO TO 30
  100 IF (L-M+1.EQ.0) THEN
        GO TO 110
      ELSE
        GO TO 105
      END IF
  105 L=L+1
      GO TO 25
  110 IF (IND-1.EQ.0) THEN
        GO TO 115
      ELSE
        GO TO 120
      END IF
  115 IND=0
      GO TO 20
  120 IF (THR-ENUF.LE.0) THEN
        GO TO 200
      ELSE
        GO TO 15
      END IF
  200 IF(IC.EQ.0)GO TO 230
      DO 225 I=1,N
      DO 226 J=I,N
      IF (EIG(I)-EIG(J).LE.0) THEN
        GO TO 226
      ELSE
        GO TO 210
      END IF
  210 X=EIG(I)
      EIG(I)=EIG(J)
      EIG(J)=X
      DO 220 K=1,N
      Y=B(K,I)
      B(K,I)=B(K,J)
      B(K,J)=Y
  220 CONTINUE
  226 CONTINUE
  225 CONTINUE
  230 CONTINUE
      RETURN
      END
