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
      SUBROUTINE FIJTD(INTSYM,INDX,C1,C2,TDMO)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),INDX(*),
     *          C1(*),C2(*),TDMO(NBAST,NBAST)
*
      JSYM(L)=JSUNP(INTSYM,L)
*------
* POW: Unnecessary but warning stopping initializations
      ni=-1234567
      nk=-1234567
*------
      ICHK=0
      IADD10=IAD10(8)
100   CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 100
      IF(LEN.LT.0) RETURN
      DO 10 IN=1,LEN
      IND=ICOP1(IN)
      IF(ICHK.NE.0) THEN
        ICHK=0
        INDI=IND
*        NI=MOD(INDI,2**10)
*        NK=MOD(INDI/2**10,2**10)
        NI=IBITS(INDI, 0,10)
        NK=IBITS(INDI,10,10)
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
      IF(IVL.NE.IVVER)GO TO 13
      TDMO(NI,NK)=TDMO(NI,NK)+COP(IN)*C1(IC1)*C2(IC2)
      IF(NI.NE.NK) TDMO(NK,NI)=TDMO(NK,NI)+COP(IN)*C2(IC1)*C1(IC2)
      GO TO 10
13    INDA=IRC(IVL)+IC1
      INDB=IRC(IVL)+IC2
      NA=INDX(INDA)
      NB=INDX(INDB)
      NS1=JSYM(INDA)
      NS1L=MUL(NS1,LSYM)
      INUM=NVIR(NS1L)
      IF(IVL.GE.2)INUM=NVPAIR(NS1L)
      TERM=DDOT_(INUM,C1(NA+1),1,C2(NB+1),1)
      TDMO(NI,NK)=TDMO(NI,NK)+COP(IN)*TERM
      IF(NI.EQ.NK) GOTO 10
      TERM=DDOT_(INUM,C2(NA+1),1,C1(NB+1),1)
      TDMO(NK,NI)=TDMO(NK,NI)+COP(IN)*TERM
10    CONTINUE
      GO TO 100
      END
