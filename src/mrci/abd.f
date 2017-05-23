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
      SUBROUTINE ABD(ICSPCK,INTSYM,INDX,C,DMO,A,B,F,JREFX)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION ICSPCK(*),INTSYM(*),INDX(*),
     *          C(*),DMO(*),A(*),B(*),
     *          F(*),JREFX(*)
      DIMENSION IPOA(9),IPOF(9)
      DIMENSION IOC(55)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CPAM97      JO(L)=UNPACK(CSPCK((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JO(L)=ICUNP(ICSPCK,L)
CPAM96      JSYM(L)=UNPACK(INTSYM((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(INTSYM,L)
C SCRATCH SPACE: A(),B(),F().
      CALL CSCALE(INDX,INTSYM,C,SQ2)
      NCLIM=4
      IF(IFIRST.NE.0)NCLIM=2
      ENPINV=1.0D00/ENP
C MOVE DENSITY MATRIX TO F IN SYMMETRY BLOCKS
      CALL IPO(IPOF,NVIR,MUL,NSYM,1,-1)
      DO 10 IASYM=1,NSYM
        IAB=IPOF(IASYM)
        NA1=NVIRP(IASYM)+1
        NA2=NVIRP(IASYM)+NVIR(IASYM)
        DO 15 NA=NA1,NA2
          DO 20 NB=NA1,NA2
            IAB=IAB+1
            NAB=IROW(LN+NA)+LN+NB
            IF(NB.GT.NA)NAB=IROW(LN+NB)+LN+NA
            F(IAB)=DMO(NAB)
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
        TSUM=ENPINV*C(INDA)**2
        GO TO 106
      END IF
      MYSYM=JSYM(INDA)
      MYL=MUL(MYSYM,LSYM)
      INMY=INDX(INDA)+1
      IF(INDA.GT.IRC(2))GO TO 25
C DOUBLET-DOUBLET INTERACTIONS
      IF(NVIR(MYL).EQ.0)GO TO 40
      CALL FMUL2(C(INMY),C(INMY),A,NVIR(MYL),NVIR(MYL),1)
      IPF=IPOF(MYL)+1
      IN=IPOF(MYL+1)-IPOF(MYL)
      CALL DAXPY_(IN,ENPINV,A,1,F(IPF),1)
      NVIRA=NVIR(MYL)
      LNA=LN+NVIRP(MYL)
      IIA=IROW(LNA+1)
      TSUM=0.0D00
      DO 130 I=1,NVIRA
        SUM=ENPINV*C(INMY)**2
        INMY=INMY+1
        TSUM=TSUM+SUM
        IIA=IIA+LNA+I
        DMO(IIA)=DMO(IIA)+SUM
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
      IF(NVIR(ICSYM).EQ.0)GO TO 70
      IF(MYL.NE.1) THEN
        IF(IASYM.GT.ICSYM) THEN
          CALL MTRANS(C(INMY+IPOA(IASYM)),1,A,1,NVIR(IASYM),NVIR(ICSYM))
        ELSE
          NAC=NVIR(IASYM)*NVIR(ICSYM)
          IF(IFT.EQ.0)CALL DCOPY_(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
          IF(IFT.EQ.1)CALL VNEG(C(INMY+IPOA(ICSYM)),1,A,1,NAC)
        END IF
      ELSE
        IF(IFT.EQ.0)CALL SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
        IF(IFT.EQ.1)CALL SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
      END IF
      NVIRA=NVIR(IASYM)
      NVIRC=NVIR(ICSYM)
      CALL FMUL2(A,A,B,NVIR(IASYM),NVIR(IASYM),NVIR(ICSYM))
      IPF=IPOF(IASYM)+1
      CALL DAXPY_(IAB,ENPINV,B,1,F(IPF),1)
      INN=1
      LNC=LN+NVIRP(ICSYM)
      IIC=IROW(LNC+1)
      DO 105 I=1,NVIRC
        SUM=ENPINV*DDOT_(NVIRA,A(INN),1,A(INN),1)
        TSUM=TSUM+SUM
        IIC=IIC+LNC+I
        DMO(IIC)=DMO(IIC)+SUM
        INN=INN+NVIRA
105   CONTINUE
70    CONTINUE
      TSUM=TSUM/2
106   CONTINUE
      IJ=0
      DO 107 I=1,LN
        IJ=IJ+I
        DMO(IJ)=DMO(IJ)+IOC(I)*TSUM
107   CONTINUE
40    CONTINUE
      DO 410 IASYM=1,NSYM
        IAB=IPOF(IASYM)
        NA1=NVIRP(IASYM)+1
        NA2=NVIRP(IASYM)+NVIR(IASYM)
        DO 415 NA=NA1,NA2
          DO 420 NB=NA1,NA2
            IAB=IAB+1
            IF(NA.GE.NB) GOTO 420
            NAB=IROW(LN+NB)+LN+NA
            DMO(NAB)=F(IAB)
420       CONTINUE
415     CONTINUE
410   CONTINUE
      TR=0.0D00
      IJ=0
      DO 510 I=1,NCSHT
        IJ=IJ+I
        TR=TR+DMO(IJ)
510   CONTINUE
      IF(ABS(TR-DBLE(NELEC)) .GT. 1.0D-8) WRITE(6,310)TR
310   FORMAT(/,6X,'TRACE OF DENSITY MATRIX',F16.8)
      CALL CSCALE(INDX,INTSYM,C,SQ2INV)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(JREFX)
      END
