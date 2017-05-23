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
      SUBROUTINE AB_MRCI(ICSPCK,INTSYM,INDX,C,S,FC,A,B,FK)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "mrci.fh"
      DIMENSION ICSPCK(*),INTSYM(*),INDX(*),
     *          C(*),S(*),FC(*),A(*),B(*),
     *          FK(*)
      DIMENSION IPOA(9),IPOF(9)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CPAM96      JSYM(L)=UNPACK(INTSYM((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(INTSYM,L)
      CALL CSCALE(INDX,INTSYM,C,SQ2)
      CALL CSCALE(INDX,INTSYM,S,SQ2INV)
      NCLIM=4
      IF(IFIRST.NE.0)NCLIM=2
C MOVE FOCK MATRIX TO FK IN SYMMETRY BLOCKS
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
            FK(IAB)=FC(NAB)
            IF(NA.EQ.NB)FK(IAB)=0.0D00
20        CONTINUE
15      CONTINUE
10    CONTINUE
      ITAIL=IRC(NCLIM)
      DO 40 INDA=1,ITAIL
      IF(INDA.LE.IRC(1))GO TO 40
      MYSYM=JSYM(INDA)
      MYL=MUL(MYSYM,LSYM)
      INMY=INDX(INDA)+1
      IF(INDA.GT.IRC(2))GO TO 25
C DOUBLET-DOUBLET INTERACTIONS
      IF(NVIR(MYL).NE.0) THEN
        CALL VCLR(A,1,NVIR(MYL))
        CALL FMMM(FK(IPOF(MYL)+1),C(INMY),A,NVIR(MYL),1,NVIR(MYL))
        CALL DAXPY_(NVIR(MYL),1.0D00,A,1,S(INMY),1)
      END IF
      GO TO 40
C TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
25    IFT=1
      IF(INDA.GT.IRC(3))IFT=0
      CALL IPO(IPOA,NVIR,MUL,NSYM,MYL,IFT)
CPAM97      IN=0
CPAM97      TSUM=0.0D00
      DO 70 IASYM=1,NSYM
      IAB=IPOF(IASYM+1)-IPOF(IASYM)
      IF(IAB.EQ.0)GO TO 70
      ICSYM=MUL(MYL,IASYM)
      IF(NVIR(ICSYM).EQ.0)GO TO 70
      IF(MYL.NE.1)GO TO 30
      IF(IFT.EQ.0)CALL SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
      IF(IFT.EQ.1)CALL SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
      NAA=NVIR(IASYM)*NVIR(IASYM)
      CALL VCLR(B,1,NAA)
      CALL FMMM(FK(IPOF(IASYM)+1),A,B,NVIR(IASYM),NVIR(IASYM),
     *NVIR(IASYM))
      CALL VCLR(A,1,NAA)
      CALL DAXPY_(NAA,1.0D00,B,1,A,1)
      IF(IFT.EQ.1)GO TO 230
      CALL SIADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
      CALL VCLR(A,1,NAA)
      GO TO 70
230   CALL TRADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
      CALL VCLR(A,1,NAA)
      GO TO 70
30    NAC=NVIR(IASYM)*NVIR(ICSYM)
      CALL VCLR(A,1,NAC)
      IF(IASYM.GT.ICSYM)GO TO 31
      CALL FMMM(FK(IPOF(IASYM)+1),C(INMY+IPOA(ICSYM)),A,
     *NVIR(IASYM),NVIR(ICSYM),NVIR(IASYM))
      CALL DAXPY_(NAC,1.0D00,A,1,S(INMY+IPOA(ICSYM)),1)
      GO TO 70
31    CALL FMMM(C(INMY+IPOA(IASYM)),FK(IPOF(IASYM)+1),A,
     *NVIR(ICSYM),NVIR(IASYM),NVIR(IASYM))
      CALL DAXPY_(NAC,1.0D00,A,1,S(INMY+IPOA(IASYM)),1)
      GO TO 70
70    CONTINUE
40    CONTINUE
      CALL CSCALE(INDX,INTSYM,C,SQ2INV)
      CALL CSCALE(INDX,INTSYM,S,SQ2)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(ICSPCK)
      END
