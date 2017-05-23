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
      SUBROUTINE AT34R(N,ISIZE,CHARGE,
     *                 SMAT,V,H,EV2,MULT,BU,P,G,EIG,SINV,REVT,
     *                 AUX,OVE,EW,E,AA,RR,TT,iprint,
     *                 VEXTT,PVPT,EVN1,EVN2,RE1R,AUXI,W1W1,W1E0W1)
C
C     INPUT: SMAT    OVERLAP MATRIX
C            V     POTENTIAL
C            H     RELATIVISTIC KINETIC ENERGY
C            EV2   PVP INTEGRALS
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "RelLight.fh"
      DIMENSION V(ISIZE),SMAT(ISIZE),MULT(ISIZE),P(ISIZE),G(ISIZE),
     *          H(ISIZE),BU(ISIZE),EV2(ISIZE),
     *          EIG(N,N),SINV(N,N),REVT(N,N),AUX(N,N),OVE(N,N),
     *          EW(N),E(N),AA(N),RR(N),TT(N)
#include "relmp.fh"
      DIMENSION VEXTT(ISIZE),PVPT(ISIZE)
      DIMENSION EVN1(N,N),EVN2(N,N)
      DIMENSION RE1R(N,N)
      DIMENSION AUXI(N,N)
      DIMENSION W1W1(N,N),W1E0W1(N,N)
*
C      CALL PRMAT(6,SMAT,N,0,'SMAT    ')
      VELIT=CLight
      ISIZE=N*(N+1)/2
      TOL=1.D-14
      PREA=1/(VELIT*VELIT)
      CON2=PREA+PREA
      CON=1.D0/PREA
      MULT(1)=0
      DO 20 I=1,N
20    MULT(I+1)=MULT(I)+I
C
C     SCHMIDT-ORTHOGONALIZE OVERLAP MATRIX
C
      CALL SOG(N,SMAT,SINV,P,OVE,EW)
      CALL FILLMA(N,SMAT,OVE)
C
C--------------------------------------------------------------------
C     MATRIX REPRESENTATION CALCULATED FROM NONRELATIVISTIC T MATRIX
C--------------------------------------------------------------------
      CALL DIAG_DKH(H,N,EIG,EW,SINV,AUX,0)
      if(iprint.ge.10) then
         write(6,*) ' eigenvalues in at34r'
         write(6,*) (ew(i),i=1,n)
      endif
      DO 4 I=1,N
C
C     IF T SUFFICIENTLY SMALL, USE SERIES EXPANSION TO AVOID CANCELLATIO
C
      TT(I)=EW(I)
      RATIO=EW(I)/VELIT
      IF (RATIO.GT.0.02D0) GOTO 11
      TV1=EW(I)
      TV2=-TV1*EW(I)*PREA/2.D0
      TV3=-TV2*EW(I)*PREA
      TV4=-TV3*EW(I)*PREA*1.25D0
      EW(I)=TV1+TV2+TV3+TV4
      GOTO 12
11    EW(I)=CON*(sqrt(1.D0+CON2*EW(I))-1.D0)
12    CONTINUE
4     E(I)=EW(I)+CON
C---------------------------------------------------------------------
C     CALCULATE REVERSE TRANSFORMATION
C---------------------------------------------------------------------
C
C     CALCULATE TRANSFORMATION MATRICES
C
      DO 3 I=1,N
      DO 3 J=1,N
      AUX(I,J)=0.D0
      DO 3 K=I,N
    3 AUX(I,J)=AUX(I,J)+SINV(I,K)*EIG(K,J)
      DO 6 I=1,N
      DO 6 J=1,N
      REVT(I,J)=0.D0
      DO 5 K=1,N
    5 REVT(I,J)=REVT(I,J)+OVE(I,K)*AUX(K,J)
    6 CONTINUE
      IJ=0
      DO 8 I=1,N
      DO 8 J=1,I
      IJ=IJ+1
      H(IJ)=0.D0
      DO 7 K=1,N
    7 H(IJ)=H(IJ)+REVT(I,K)*REVT(J,K)*EW(K)
    8 CONTINUE
      IF(IRELMP.NE.11) THEN
      DO 362 I=1,N
      AA(I)=sqrt((CON+E(I)) / (2.D0*E(I)))
362   RR(I)=sqrt(CON)/(CON+E(I))
      ELSE IF(IRELMP.EQ.11) THEN                          ! RESC
      DO I=1,N
        AA(I)=(sqrt(1.0D0+CON*TT(I)*2.0D0/((CON+E(I))*(CON+E(I)))))
     &                              /(CON+E(I))           ! O OPERATOR
        RR(I)=sqrt(CON)/(CON+E(I))                       ! Q OPERATOR
      ENDDO
      ENDIF
C
C     BEYOND THIS POINT, OVE IS USED AS SCRATCH ARRAY
C
C
C    TRANSFORM V TO T-BASIS
C
      CALL TRSM_DKH(V,SINV,G,N,AUX,OVE)
      CALL TRSM_DKH(G,EIG,BU,N,AUX,OVE)
C
C    MULTIPLY
C
      IF(IRELMP.NE.11) THEN
*
      IJ=0
      DO 2005 I=1,N
      DO 2006 J=1,I
      IJ=IJ+1
      P(IJ)=-BU(IJ)*CHARGE
      VEXTT(IJ)=P(IJ) ! KEEP T-BASIS VEXT INTO VEXTT FOR HIGHER-ORDER DK
      BU(IJ)= P(IJ)*AA(I)*AA(J)
      EVN1(I,J)=BU(IJ)
      EVN1(J,I)=EVN1(I,J)
2006  CONTINUE
2005  CONTINUE
*
      ELSE IF(IRELMP.EQ.11) THEN
*
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          P(IJ)=-BU(IJ)*CHARGE
          BU(IJ)=VELIT* P(IJ)*(sqrt(RR(I)*RR(J))*AA(I)/AA(J)
     $                        +sqrt(RR(J)*RR(I))*AA(J)/AA(I))
        ENDDO
      ENDDO
*
      ENDIF
*
      CALL TRSMT(BU,REVT,V,N,AUX,OVE)
C
C     PVP INTEGRALS
C
      CALL TRSM_DKH(EV2,SINV,G,N,AUX,OVE)
      CALL TRSM_DKH(G,EIG,BU,N,AUX,OVE)
C
C    MULTIPLY
C
      IF(IRELMP.NE.11) THEN
      IJ=0
      DO 3005 I=1,N
      DO 3006 J=1,I
      IJ=IJ+1
      G(IJ)=-BU(IJ)*CHARGE
      PVPT(IJ)=G(IJ)   ! KEEP T-BASIS PVP INTO PVPT FOR HIGHER-ORDER DK
      BU(IJ)= G(IJ)*AA(I)*RR(I)*AA(J)*RR(J)
      EVN1(I,J)=EVN1(I,J)+BU(IJ)
      EVN1(J,I)=EVN1(I,J)
3006  CONTINUE
3005  CONTINUE
      ELSE IF(IRELMP.EQ.11) THEN
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          G(IJ)=-BU(IJ)*CHARGE
          BU(IJ)= G(IJ)*(RR(I)*RR(J)*AA(I)/AA(J)
     $                  +RR(J)*RR(I)*AA(J)/AA(I))*0.5D0
        ENDDO
      ENDDO
      ENDIF
      CALL TRSMT(BU,REVT,EV2,N,AUX,OVE)
C@    CALL PRMAT(6,EV2,N,0,'PVPFULL ')
      CALL ADDMA(ISIZE,EV2,V)
*
      IF(IRELMP.EQ.1.OR.IRELMP.EQ.11) GOTO 1000
C
C     CALCULATE EVEN2 OPERATOR
C
C@    CALL PRMAT(6,TT,N,1,'TT      ')
C@    CALL PRMAT(6,E,N,1,'E       ')
      CALL EVEN2(N,P,G,E,AA,RR,TT,EIG,AUX,OVE,
     &           W1W1,W1E0W1,EVN2)
*
*    TRANSFORM BACK
*
      CALL TRSMT(G,REVT,EV2,N,AUX,OVE)
      CALL ADDMA(ISIZE,EV2,V)
*
      IF(IRELMP.EQ.0.OR.IRELMP.EQ.2) GOTO 1000   ! DK2
*
C
C     ----- CALCULATE Even3r OPERATOR -----
C
      CALL EVEN3(N,P,G,E,AA,RR,TT,EIG,AUX,OVE,
     &           EVN1,VEXTT,PVPT,RE1R,W1W1,AUXI)
C
C     ------- TRANSFORM BACK FOR DK3
C
      CALL TRSMT(G,REVT,EV2,N,AUX,OVE)
      CALL ADDMA(ISIZE,EV2,V)
      IF(IRELMP.EQ.3) GOTO 1000   ! DK3
*
*     More to come here
*
 1000 CONTINUE
*
      CR=1/CHARGE
      DO 940 I=1,ISIZE
940   V(I)=-V(I)*CR
C@    CALL PRMAT(6,BU,N,0,'TOTAL H ')
      RETURN
      END
