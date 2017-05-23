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
      SUBROUTINE FIJD(INTSYM,INDX,C,DMO,JREFX,AREF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),INDX(*),
     *          C(*),DMO(*),JREFX(*),AREF(*)
*
      JSYM(L)=JSUNP(INTSYM,L)
*
      ICHK=0
      IK=0
      ENPINV=1.0D00/ENP
      DFACR=1.0D00-ENPINV
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
      IF(IVL.NE.IVVER)GO TO 13
      DMO(IK)=DMO(IK)+COP(IN)*C(IC1)*C(IC2)*ENPINV
      IF(ICPF.EQ.0)GO TO 10
      IRC1=JREFX(IC1)
      IF(IRC1.EQ.0)GO TO 10
      IRC2=JREFX(IC2)
      IF(IRC2.EQ.0)GO TO 10
      DMO(IK)=DMO(IK)+COP(IN)*AREF(IRC1)*AREF(IRC2)*(1.0D00-ENPINV)
      GO TO 10
13    INDA=IRC(IVL)+IC1
      INDB=IRC(IVL)+IC2
      NA=INDX(INDA)
      NB=INDX(INDB)
      NS1=JSYM(INDA)
      NS1L=MUL(NS1,LSYM)
      INUM=NVIR(NS1L)
      IF(IVL.GE.2)INUM=NVPAIR(NS1L)
      TERM=DDOT_(INUM,C(NA+1),1,C(NB+1),1)
      DMO(IK)=DMO(IK)+COP(IN)*TERM*ENPINV
10    CONTINUE
      GO TO 100
      END
