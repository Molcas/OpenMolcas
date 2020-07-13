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
* Copyright (C) 1986, Bernd Artur Hess                                 *
************************************************************************
      SUBROUTINE EVEN2(N,V,G,E,A,R,TT,AUXF,AUXG,AUXH,
     &                 W1W1,W1E0W1,EVN2)
C
C     EVEN2 - BERND HESS - V 1.0 - 5.2.86
C     CALCULATE EVEN2 OPERATORS
C
C
C     N       DIMENSION OF MATRICES
C     V       POTENTIAL MATRIX
C     G       MATRIX OF PVP OPERATOR. WILL CONTAIN EVEN2 OPERATORS
C             ON OUTPUT
C     E       RELATIVISTIC ENERGY (DIAGONAL)
C     A       A-FACTORS (DIAGONAL)
C     R       R-FACTORS (DIAGONAL)
C     TT      NONREL. KINETIC ENERGY (DIAGONAL)
C     AUXF,AUXG,AUXH  SCRATCH ARAYS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(N*(N+1)/2),G(N*(N+1)/2),E(N),R(N),A(N),TT(N),
     &          AUXF(N,N),AUXG(N,N),AUXH(N,N)
      DIMENSION W1W1(N,N),W1E0W1(N,N)
      DIMENSION EVN2(N,N)
*
      M=N
      IJ=0
      DO 600 I=1,N
      DO 601 J=1,I
      IJ=IJ+1
      AUXH(I,J)=0.D0
      AUXH(J,I)=0.D0
      V(IJ)=V(IJ)/(E(I)+E(J))
      G(IJ)=G(IJ)/(E(I)+E(J))
      AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)*A(J)
      AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)*A(I)
      AUXG(I,J)=R(I)*V(IJ)*A(J)
      AUXG(J,I)=R(J)*V(IJ)*A(I)
601   CONTINUE
600   CONTINUE
C
C     ARQA ARQA
C
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IF (IE.NE.0) Call SysHalt('relint')
CCC   CALL PRSQ('AUXH   1',AUXH,N)
      IJ=0
      DO 602 I=1,N
      DO 603 J=1,I
      IJ=IJ+1
      AUXG(I,J)=-0.5D0/TT(I)*G(IJ)*A(J)*R(J)
      AUXG(J,I)=-0.5D0/TT(J)*G(IJ)*A(I)*R(I)
603   CONTINUE
602   CONTINUE
C
C     ARQA AQRA
C
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
CCC   CALL PRSQ('AUXH   2',AUXH,N)
      IJ=0
      DO 604 I=1,N
      DO 605 J=1,I
      IJ=IJ+1
      AUXF(I,J)=A(I)*V(IJ)*A(J)*A(J)*R(J)
      AUXF(J,I)=A(J)*V(IJ)*A(I)*A(I)*R(I)
      AUXG(I,J)=-2.D0*TT(I)*R(I)*V(IJ)*A(J)
      AUXG(J,I)=-2.D0*TT(J)*R(J)*V(IJ)*A(I)
605   CONTINUE
604   CONTINUE
C
C     AQRA ARQA
C
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
CCC   CALL PRSQ('AUXH   3',AUXH,N)
      IJ=0
      DO 606 I=1,N
      DO 607 J=1,I
      IJ=IJ+1
      AUXG(I,J)=G(IJ)*A(J)*R(J)
      AUXG(J,I)=G(IJ)*A(I)*R(I)
607   CONTINUE
606   CONTINUE
C
C     AQRA AQRA
C
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
C
C     KEEP W1*W1 FOR HIGHER-ORDER DK
C
      DO I=1,N
        DO J=1,N
          W1W1(I,J)=AUXH(I,J)
        ENDDO
      ENDDO
C
C     1/2 EW*W + 1/2 W*WE
C
      DO 610 I=1,N
      DO 611 J=1,N
      AUXH(I,J)=0.5D0*( AUXH(I,J)*E(I) + AUXH(I,J)*E(J) )
611   CONTINUE
610   CONTINUE
C
      IJ=0
      DO 400 I=1,N
      DO 401 J=1,I
      IJ=IJ+1
      AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)*E(J)*A(J)
      AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)*E(I)*A(I)
      AUXG(I,J)=R(I)*V(IJ)*A(J)
      AUXG(J,I)=R(J)*V(IJ)*A(I)
401   CONTINUE
400   CONTINUE
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IJ=0
      DO 402 I=1,N
      DO 403 J=1,I
      IJ=IJ+1
      AUXG(I,J)=-0.5D0/TT(I)*G(IJ)*A(J)*R(J)
      AUXG(J,I)=-0.5D0/TT(J)*G(IJ)*A(I)*R(I)
403   CONTINUE
402   CONTINUE
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
CCC   CALL PRSQ('AUXH   6',AUXH,N)
      IJ=0
      DO 404 I=1,N
      DO 405 J=1,I
      IJ=IJ+1
      AUXF(I,J)=A(I)*V(IJ)*R(J)*A(J)*E(J)*A(J)
      AUXF(J,I)=A(J)*V(IJ)*R(I)*A(I)*E(I)*A(I)
      AUXG(I,J)=-2.D0*TT(I)*R(I)*V(IJ)*A(J)
      AUXG(J,I)=-2.D0*TT(J)*R(J)*V(IJ)*A(I)
405   CONTINUE
404   CONTINUE
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
      IJ=0
      DO 406 I=1,N
      DO 407 J=1,I
      IJ=IJ+1
      AUXG(I,J)=G(IJ)*A(J)*R(J)
      AUXG(J,I)=G(IJ)*A(I)*R(I)
407   CONTINUE
406   CONTINUE
      CALL CPLAB(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
CCC   CALL PRSQ('AUXH   8',AUXH,N)
C
C     SYMMETRISIEREN
C
      IJ=0
      DO 430 I=1,N
      DO 431 J=1,I
      IJ=IJ+1
      G(IJ)=-0.5D0*(AUXH(I,J)+AUXH(J,I))
 431  CONTINUE
 430  CONTINUE
CCC   CALL PRM('OUTPUT  ',G,N)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_real_array(W1E0W1)
        CALL Unused_real_array(EVN2)
      END IF
      END
