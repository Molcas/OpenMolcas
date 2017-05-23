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
      SUBROUTINE EVEN3(N,V,G,E,A,R,TT,AUXF,AUXG,AUXH,
     $                 EVEN1,VEXTT,PVPT,RE1R,W1W1,AUXI)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(N*(N+1)/2),G(N*(N+1)/2),E(N),R(N),A(N),TT(N),
     &          AUXF(N,N),AUXG(N,N),AUXH(N,N)
      DIMENSION EVEN1(N,N)
      DIMENSION RE1R(N,N)
      DIMENSION VEXTT(N*(N+1)/2),PVPT(N*(N+1)/2)
      DIMENSION AUXI(N,N)
      DIMENSION W1W1(N,N)
C
C     ----- CONSTRUCT RE1R -----
C
      CALL DKRE1R(A,R,E,TT,V,G,RE1R,VEXTT,PVPT,N)
C
      M=N
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXH(I,J)=0.D0
          AUXH(J,I)=0.D0
          V(IJ)=VEXTT(IJ)/(E(I)+E(J))
          G(IJ)=PVPT(IJ)/(E(I)+E(J))
        ENDDO
      ENDDO
C
C     ----- W1*W1*E1 -----
C
C     ------- 1/2 E1*W1*W1 + 1/2 W1*W1*E1
C
      DO I=1,N
        DO J=1,N
          AUXF(I,J)=0.5D0*EVEN1(I,J)
        ENDDO
      ENDDO
      CALL CpLabr(W1W1,AUXF,N,N,N,M,M,AUXH,M,IE)
      CALL CpLabr(AUXF,W1W1,N,N,N,M,M,AUXH,M,IE)
C
C     ----- W1*E1*W1 TERM -----
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)/R(J)/TT(J)*0.5D0
          AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)/R(I)/TT(I)*0.5D0
          AUXG(I,J)=-A(I)*V(IJ)*A(J)
          AUXG(J,I)=-A(J)*V(IJ)*A(I)
        ENDDO
      ENDDO
      CALL CZERO2(AUXI,N,N,N)
      CALL CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
      CALL CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXF(I,J)=A(I)*V(IJ)*A(J)
          AUXF(J,I)=A(J)*V(IJ)*A(I)
          AUXG(I,J)=-A(I)/R(I)*G(IJ)*A(J)*R(J)/TT(I)*0.5D0
          AUXG(J,I)=-A(J)/R(J)*G(IJ)*A(I)*R(I)/TT(J)*0.5D0
        ENDDO
      ENDDO
      CALL CZERO2(AUXI,N,N,N)
      CALL CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
      CALL CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXF(I,J)=A(I)*V(IJ)*A(J)
          AUXF(J,I)=A(J)*V(IJ)*A(I)
          AUXG(I,J)=A(I)*V(IJ)*A(J)
          AUXG(J,I)=A(J)*V(IJ)*A(I)
        ENDDO
      ENDDO
      CALL CZERO2(AUXI,N,N,N)
      CALL CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
      CALL CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)/R(J)/TT(J)*0.5D0
          AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)/R(I)/TT(I)*0.5D0
          AUXG(I,J)=A(I)/R(I)*G(IJ)*A(J)*R(J)/TT(I)*0.5D0
          AUXG(J,I)=A(J)/R(J)*G(IJ)*A(I)*R(I)/TT(J)*0.5D0
        ENDDO
      ENDDO
      CALL CZERO2(AUXI,N,N,N)
      CALL CpLabr(AUXF,RE1R,N,N,N,M,M,AUXI,M,IE)
      CALL CpLabr(AUXI,AUXG,N,N,N,M,M,AUXH,M,IE)
C
C     ----- SYMMETRISIEREN -----
C
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          G(IJ)=0.5D0*(AUXH(I,J)+AUXH(J,I))
        ENDDO
      ENDDO
C
      RETURN
      END
