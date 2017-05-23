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
      SUBROUTINE FIJ_MRCI(ICSPCK,INTSYM,INDX,C,S,FC,A,B,FK,DBK)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION ICSPCK(*),INTSYM(*),INDX(*),
     *          C(*),S(*),FC(*),A(*),B(*),
     *          FK(*),DBK(*)
*
      JSYM(L)=JSUNP(INTSYM,L)
*
      ICHK=0
      IK=0
      IADD25=0
      CALL dDAFILE(Lu_25,2,FC,NBTRI,IADD25)
      IADD10=IAD10(8)
100   CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 100
      IF(LEN.LT.0)GO TO 200
      DO 10 IN=1,LEN
      IND=ICOP1(IN)
      IF(ICHK.NE.0) THEN
        ICHK=0
        INDI=IND
*        NI=MOD(INDI,2**10)
*        NK=MOD(INDI/2**10,2**10)
        NI=IBITS(INDI, 0,10)
        NK=IBITS(INDI,10,10)
        IK=IROW(NK)+NI
        GO TO 10
      END IF
      IF(IND.EQ.0) THEN
        ICHK=1
        GO TO 10
      END IF
*      IVL=MOD(IND,2**6)
*      IC2=MOD(IND/2**6,2**13)
*      IC1=MOD(IND/2**19,2**13)
      IVL=IBITS(IND, 0, 6)
      IC2=IBITS(IND, 6,13)
      IC1=IBITS(IND,19,13)
      COPI=COP(IN)*FC(IK)
      IF(IVL.EQ.IVVER) THEN
        S(IC1)=S(IC1)+COPI*C(IC2)
        S(IC2)=S(IC2)+COPI*C(IC1)
        GO TO 10
      END IF
      INDA=IRC(IVL)+IC1
      INDB=IRC(IVL)+IC2
      NA=INDX(INDA)
      NB=INDX(INDB)
      NS1=JSYM(INDA)
      NS1L=MUL(NS1,LSYM)
      INUM=NVIR(NS1L)
      IF(IVL.GE.2)INUM=NVPAIR(NS1L)
      CALL DAXPY_(INUM,COPI,C(NB+1),1,S(NA+1),1)
      CALL DAXPY_(INUM,COPI,C(NA+1),1,S(NB+1),1)
10    CONTINUE
      GO TO 100
200   IF(ITER.EQ.0)RETURN
      CALL AI_MRCI(INTSYM,INDX,C,S,FC,C,C,A,B,FK,DBK,0)
      IF(ITER.EQ.1 .AND. IREST.EQ.0)RETURN
      CALL AB_MRCI(ICSPCK,INTSYM,INDX,C,S,FC,A,B,FK)
      RETURN
      END
