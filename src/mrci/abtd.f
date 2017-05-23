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
      SUBROUTINE ABTD(ICSPCK,INTSYM,INDX,C1,C2,TDMO,A1,A2,F)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION ICSPCK(*),INTSYM(*),INDX(*),
     *          C1(*),C2(*),TDMO(NBAST,NBAST),
     *          A1(*),A2(*),F(*)
      DIMENSION IPOA(9),IPOF(9)
      DIMENSION IOC(55)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CPAM97      JO(L)=UNPACK(CSPCK((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JO(L)=ICUNP(ICSPCK,L)
CPAM96      JSYM(L)=UNPACK(INTSYM((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(INTSYM,L)
C CALCULATE A) TRANSITION DENSITY ELEMENTS OF TYPE TDMO(A,B)
C           B) DIAGONAL ELEMENTS TDMO(I,I) AND TDMO(A,A)
C SCRATCH SPACES: A1(),A2(), SIZE NEEDED IS NVMAX**2
C                  ALSO F(), SIZE NEEDED IS NVSQ
      CALL CSCALE(INDX,INTSYM,C1,SQ2)
      CALL CSCALE(INDX,INTSYM,C2,SQ2)
      NCLIM=4
      IF(IFIRST.NE.0)NCLIM=2
C MOVE TRANSITION DENSITY MATRIX TO F IN SYMMETRY BLOCKS
      CALL IPO(IPOF,NVIR,MUL,NSYM,1,-1)
      DO 10 IASYM=1,NSYM
        IAB=IPOF(IASYM)
        NA1=NVIRP(IASYM)+1
        NA2=NVIRP(IASYM)+NVIR(IASYM)
        DO 15 NA=NA1,NA2
          DO 20 NB=NA1,NA2
            IAB=IAB+1
            F(IAB)=TDMO(LN+NA,LN+NB)
20        CONTINUE
15      CONTINUE
10    CONTINUE
      II1=0
      ITAIL=IRC(NCLIM)
      DO 40 INDA=1,ITAIL
        DO 110 I=1,LN
          II1=II1+1
          IOC(I)=(1+JO(II1))/2
110     CONTINUE
      IF(INDA.LE.IRC(1)) THEN
        TSUM=C1(INDA)*C2(INDA)
        GO TO 106
      END IF
      MYSYM=JSYM(INDA)
      MYL=MUL(MYSYM,LSYM)
      INMY=INDX(INDA)+1
      IF(INDA.GT.IRC(2))GO TO 25
C DOUBLET-DOUBLET INTERACTIONS
      IF(NVIR(MYL).EQ.0)GO TO 40
      IPF=IPOF(MYL)+1
      IN=IPOF(MYL+1)-IPOF(MYL)
      NVIRA=NVIR(MYL)
      CALL DGER(NVIRA,NVIRA,1.0D00,C1(INMY),1,
     *           C2(INMY),1,F(IPF),NVIRA)
      LNA=LN+NVIRP(MYL)
      TSUM=0.0D00
      DO 130 I=1,NVIRA
        TERM=C1(INMY-1+I)*C2(INMY-1+I)
        IA=LNA+I
        TDMO(IA,IA)=TDMO(IA,IA)+TERM
        TSUM=TSUM+TERM
130   CONTINUE
      GO TO 106
C TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
25    IFT=1
      IF(INDA.GT.IRC(3))IFT=0
      CALL IPO(IPOA,NVIR,MUL,NSYM,MYL,IFT)
      IN=0
      TSUM=0.0D00
      DO 70 IASYM=1,NSYM
      IAB=IPOF(IASYM+1)-IPOF(IASYM)
      IF(IAB.EQ.0)GO TO 70
      ICSYM=MUL(MYL,IASYM)
      NVIRA=NVIR(IASYM)
      NVIRC=NVIR(ICSYM)
      IF(NVIRC.EQ.0)GO TO 70
      IF(MYL.NE.1) THEN
        IF(IASYM.GT.ICSYM) THEN
          CALL MTRANS(C1(INMY+IPOA(IASYM)),1,A1,1,NVIRA,NVIRC)
          CALL MTRANS(C2(INMY+IPOA(IASYM)),1,A2,1,NVIRA,NVIRC)
        ELSE
          NAC=NVIRA*NVIRC
          IF(IFT.EQ.0) THEN
            CALL DCOPY_(NAC,C1(INMY+IPOA(ICSYM)),1,A1,1)
            CALL DCOPY_(NAC,C2(INMY+IPOA(ICSYM)),1,A2,1)
          ELSE
            CALL VNEG(C1(INMY+IPOA(ICSYM)),1,A1,1,NAC)
            CALL VNEG(C2(INMY+IPOA(ICSYM)),1,A2,1,NAC)
          END IF
        END IF
      ELSE
        IF(IFT.EQ.0)THEN
          CALL SQUAR(C1(INMY+IPOA(IASYM)),A1,NVIRA)
          CALL SQUAR(C2(INMY+IPOA(IASYM)),A2,NVIRA)
        ELSE
          CALL SQUARM(C1(INMY+IPOA(IASYM)),A1,NVIRA)
          CALL SQUARM(C2(INMY+IPOA(IASYM)),A2,NVIRA)
        END IF
      END IF
      IPF=IPOF(IASYM)+1
      CALL DGEMM_('N','T',NVIRA,NVIRA,NVIRC,1.0D00,A1,NVIRA,
     *            A2,NVIRA,1.0D00,F(IPF),NVIRA)
      INN=1
      LNC=LN+NVIRP(ICSYM)
      DO 105 I=1,NVIRC
        TERM=DDOT_(NVIRA,A1(INN),1,A2(INN),1)
        TSUM=TSUM+TERM
        IC=LNC+I
        TDMO(IC,IC)=TDMO(IC,IC)+TERM
        INN=INN+NVIRA
105   CONTINUE
70    CONTINUE
      TSUM=TSUM/2
106   CONTINUE
      DO 107 I=1,LN
        TDMO(I,I)=TDMO(I,I)+IOC(I)*TSUM
107   CONTINUE
40    CONTINUE
      DO 410 IASYM=1,NSYM
        IAB=IPOF(IASYM)
        NA1=NVIRP(IASYM)+1
        NA2=NVIRP(IASYM)+NVIR(IASYM)
        DO 415 NA=NA1,NA2
          DO 420 NB=NA1,NA2
            IAB=IAB+1
            IF(NA.EQ.NB) GOTO 420
            TDMO(LN+NA,LN+NB)=F(IAB)
420       CONTINUE
415     CONTINUE
410   CONTINUE
      CALL CSCALE(INDX,INTSYM,C1,SQ2INV)
      CALL CSCALE(INDX,INTSYM,C2,SQ2INV)
      RETURN
      END
