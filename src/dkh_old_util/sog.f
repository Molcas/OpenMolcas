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
      SUBROUTINE SOG(N,SS,SINV,P,G,A1)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SS(N*(N+1)/2),P(N*(N+1)/2),G(N*(N+1)/2),A1(N),SINV(N,N)
C
C     SUBROUTINE TO CALCULATE TRANSFORMATION TO SCHMIDT-
C     ORTHOGONALIZED BASIS.
C     N              DIMENSION OF MATRICES. ISIZE=N*(N+1)/2
C     SS(ISIZE)      ORIGINAL OVERLAP MATRIX (LOWER TRIANGULAR)
C                    WILL NOT BE DESTROYED
C     P (ISIZE)      OUTPUT TRANSFORMATION MATRIX
C     G (ISIZE)      SCRATCH
C     A1(N)          SCRATCH
C
      JL=0
      IQ=0
      DO 340 J=1,N
         IL=JL
         JQ=IQ
         S1KK=SS(IQ+J)
         G(IL+J)=1.D0
         IF(J.EQ.1)GO TO 341
         J1=J-1
         JL=0
         DO 342 K=1,J1
            LG=JQ
            ETOT=0.D0
            DO 343 L=1,K
               LG=LG+1
               JL=JL+1
  343       ETOT=ETOT+SS(LG)*G(JL)
            S1KK=S1KK-ETOT*ETOT
  342    A1(K)=ETOT
         IF=1
         JL=IL
         DO 344 K=1,J1
         SUM=0.D0
         JL=JL+1
         IF=IF+K-1
         IH=IF
         DO 345 L=K,J1
         IH=IH+L-1
  345    SUM=SUM+A1(L)*G(IH)
  344    G(JL)=-SUM
  341    S1KK=1.D0/sqrt(S1KK)
         JL=IL
         DO 340 K=1,J
         JL=JL+1
         IQ=IQ+1
         G(JL)=G(JL)*S1KK
  340 P(IQ)=G(JL)
*
      IJ=0
      DO 1 I=1,N
         DO 2 J=1,I
            IJ=IJ+1
            SINV(I,J)=0.D0
2        SINV(J,I)=P(IJ)
1     CONTINUE
      RETURN
      END
