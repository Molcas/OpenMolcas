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
        SUBROUTINE BNDINV(       A,      EL,       N,  DETERM,   EPSIL,
     &                       ITEST,   NSIZE)
C
C       REAL*8 MATRIX INVERSION SUBROUTINE
C       FROM "DLYTAP".
C
C*      REAL*8 E,F
C*      REAL*8 A,EL,D,SQRT,C,S,DETERP
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(NSIZE,*),EL(NSIZE,*)
        INDSNL=0
        IF(N.LT.2)GO TO 140
        ISL2=0
        K000FX=2
        IF(ISL2.EQ.0)INDSNL=2
        IF(ISL2.EQ.1)INDSNL=1
C       CALL SLITET(2,INDSNL)
C       CALL OVERFL(K000FX)
C       CALL DVCHK(K000FX)
C
C       SET EL = IDENTITY MATRIX
        DO I=1,N
         DO J=1,N
          EL(I,J)=0.0D0
         END DO
         EL(I,I)=1.0D0
        END DO
C
C       TRIANGULARIZE A, FORM EL
C
        N1=N-1
        M=2
        DO J=1,N1
         DO I=M,N
          IF(A(I,J).EQ.0.0D0)GO TO 45
          D=sqrt(A(J,J)*A(J,J)+A(I,J)*A(I,J))
          C=A(J,J)/D
          S=A(I,J)/D
          DO K=J,N
           D=C*A(J,K)+S*A(I,K)
           A(I,K)=C*A(I,K)-S*A(J,K)
           A(J,K)=D
          END DO
          DO K=1,N
           D=C*EL(J,K)+S*EL(I,K)
           EL(I,K)=C*EL(I,K)-S*EL(J,K)
           EL(J,K)=D
          END DO
 45       CONTINUE
         END DO
         M=M+1
        END DO
C       CALL OVERFL(K000FX)
C       GO TO (140,51),K000FX
C
C       CALCULATE THE DETERMINANT
        CONTINUE
        DETERP=A(1,1)
        DO I=2,N
         DETERP=DETERP*A(I,I)
        END DO
        DETERM=DETERP
C       CALL OVERFL(K000FX)
C       GO TO (140,520,520),K000FX
C
C       IS MATRIX SINGULAR
c 520    CONTINUE
        F=A(1,1)
        E=A(1,1)
        DO I=2,N
         IF(abs(F).LT.abs(A(I,I)))F=A(I,I)
         IF(abs(E).GT.abs(A(I,I)))E=A(I,I)
        END DO
        EPSILP=EPSIL
        IF(EPSILP.LE.0.0D0)EPSILP=1.0D-8
        RAT=E/F
        IF(ABS(RAT).LT.EPSILP)GO TO 130
C
C       INVERT TRIANGULAR MATRIX
        J=N
        DO J1=1,N
         A(J,J)=1.0D0/A(J,J)
         I=J-1
         DO I1=2,J
          D=0.0D0
          DO K=I+1,J
            D=D+A(I,K)*A(K,J)
          END DO
          A(I,J)=-D/A(I,I)
          I=I-1
         END DO
         J=J-1
        END DO
C       CALL OVERFL(K000FX)
C       GO TO (140,103,103),K000FX

C103    CALL DVCHK(K000FX)
C       GO TO (140,105),K000FX
C
C       PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
        M=1
        DO 120 I=1,N
        DO 118 J=1,N
        D=0.0D0
        DO 107 K=M,N
        D=D+A(I,K)*EL(K,J)
 107    CONTINUE
        EL(I,J)=D
 118    CONTINUE
        M=M+1
 120    CONTINUE
C       CALL OVERFL(K000FX)
C       GO TO (140,123,123),K000FX
C
C       RECOPY EL TO A
        DO 124 I=1,N
        DO 125 J=1,N
        A(I,J)=EL(I,J)
 125    CONTINUE
 124    CONTINUE
        ITEST=0
C126    IF(INDSNL.EQ.1)CALL SLITE(2)
 126    IF(INDSNL.EQ.1)ISL2=1
        RETURN
C
 130    ITEST=1
        GO TO 126
 140    ITEST=-1
        GO TO 126
        END
