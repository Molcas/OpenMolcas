************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
cpgi$g opt=1
      SUBROUTINE FIJ(ICASE,JSY,INDEX,C,S,FC,A,B,FK,DBK,ENP,
     *EPP)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION JSY(*),INDEX(*),C(*),S(*),FC(*),A(*),B(*),
     &             FK(*),DBK(*),ENP(*),EPP(*)
      DIMENSION ICASE(*)
      PARAMETER (IPOW6=2**6, IPOW10=2**10, IPOW19=2**19)
*
      JSYM(L)=JSUNP(JSY,L)
*
      IK = 0 ! dummy initialize
      NOB2=IROW(NORBT+1)
C      IF(IDENS.EQ.1)WRITE(6,876)(FC(I),I=1,NOB2)
C  876 FORMAT(1X,'FIJ',5F12.6)
      ICHK=0
      IF(IDENS.EQ.1)GO TO 105
      NOB2=IROW(NORBT+1)
      CALL SETZ(FC,NOB2)
      IADD25=0
      CALL dDAFILE(Lu_25,2,FC,NOB2,IADD25)
      IF(ITER.EQ.1)GO TO 200
105   IADD10=IAD10(8)
100   CALL dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      CALL iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 100
      IF(LEN.LT.0)GO TO 200
      DO 10 IN=1,LEN
      IND=ICOP1(IN)
      IF(ICHK.NE.0)GO TO 460
      IF(IND.NE.0)GO TO 11
      ICHK=1
      GO TO 10
460   ICHK=0
      INDI=IND
*      NI=MOD(INDI,1024)
*      NK=MOD(INDI/IPOW10,1024)
      NI=IBITS(INDI,0,10)
      NK=IBITS(INDI,10,10)
      IK=IROW(NK)+NI
      GO TO 10
11    CONTINUE
*      IVL=MOD(IND,64)
*      IC2=MOD(IND/IPOW6,8192)
*      IC1=MOD(IND/IPOW19,8192)
      IVL=IBITS(IND, 0, 6)
      IC2=IBITS(IND,6,13  )
      IC1=IBITS(IND,19,13  )
      COPI=COP(IN)*FC(IK)
      IF(IVL.NE.IV0)GO TO 13
      IF(IC1.NE.IREF0)GO TO 16
      IF(IDENS.EQ.1)GO TO 18
      COPI=COPI/SQRT(ENP(IC2))
      S(IC2)=S(IC2)+COPI
      IF(ITER.EQ.1)GO TO 10
      EPP(IC2)=EPP(IC2)+COPI*C(IC2)
      GO TO 10
18    FC(IK)=FC(IK)+COP(IN)*C(IC1)*C(IC2)/ENP(IC2)
      GO TO 10
16    IF(IC2.NE.IREF0)GO TO 17
      IF(IDENS.EQ.1)GO TO 19
      COPI=COPI/SQRT(ENP(IC1))
      S(IC1)=S(IC1)+COPI
      IF(ITER.EQ.1)GO TO 10
      EPP(IC1)=EPP(IC1)+COPI*C(IC1)
      GO TO 10
19    FC(IK)=FC(IK)+COP(IN)*C(IC1)*C(IC2)/ENP(IC1)
      GO TO 10
17    IF(IDENS.EQ.1)GO TO 21
      FACS=D1
      S(IC1)=S(IC1)+FACS*COPI*C(IC2)
      S(IC2)=S(IC2)+FACS*COPI*C(IC1)
      GO TO 10
21    FC(IK)=FC(IK)+COP(IN)*C(IC1)*C(IC2)/
     *(SQRT(ENP(IC1))*SQRT(ENP(IC2)))
      GO TO 10
13    INDA=IRC(IVL)+IC1
      INDB=IRC(IVL)+IC2
      NA=INDEX(INDA)
      NB=INDEX(INDB)
      NS1=JSYM(INDA)
      NS1L=MUL(NS1,LSYM)
      INUM=NVIR(NS1L)
      IF(IVL.GE.2)INUM=NNS(NS1L)
      IF(IDENS.EQ.1)GO TO 15
      FACS=D1
      CALL DAXPY_(INUM,COPI*FACS,C(NB+1),1,S(NA+1),1)
      CALL DAXPY_(INUM,COPI*FACS,C(NA+1),1,S(NB+1),1)
      GO TO 10
15    TERM=DDOT_(INUM,C(NA+1),1,C(NB+1),1)
      FC(IK)=FC(IK)+COP(IN)*TERM/(SQRT(ENP(INDA))*SQRT(ENP(INDB)))
10    CONTINUE
      GO TO 100
C  200 IF(IDENS.EQ.1)WRITE(6,876)(FC(I),I=1,NOB2)
200   CALL dAI_CPF(C)
      IF(ITER.EQ.1)RETURN
      CALL AB(ICASE,JSY,INDEX,C,S,FC,A,B,FK,ENP)
      RETURN
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE dAI_CPF(C)
      USE ISO_C_BINDING
      REAL*8, TARGET :: C(*)
      INTEGER, POINTER :: iC(:)
      CALL C_F_POINTER(C_LOC(C(1)),iC,[1])
      CALL AI_CPF(JSY,INDEX,C,S,FC,C,iC,A,B,FK,DBK,ENP,EPP,0)
      NULLIFY(iC)
      END SUBROUTINE dAI_CPF
*
      END
