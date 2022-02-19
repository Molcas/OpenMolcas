!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine ABCI(JSY,INDEX,C,S,BMN,IBMN,BIAC,BICA,BUFIN)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
dimension JSY(*), index(*), C(*), S(*), BMN(*), IBMN(*), BIAC(*), BICA(*), BUFIN(*)
parameter(IPOW6=2**6,IPOW19=2**19)
! Statement function
JSYM(L) = JSUNP_CPF(JSY,L)

INUM = IRC(4)-IRC(3)
call PSQ2(C,S,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
ICHK = 0
INSIN = KBUFF1
IAD15 = IADABCI
IADD10 = IAD10(4)
call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
IN = 2
NSAVE = ICOP1(IN)
100 NI = NSAVE
IOUT = 0
110 IN = IN+1
if (IN <= LEN) GO TO 15
call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN <= 0) GO TO 5
IN = 1
15 if (ICHK /= 0) GO TO 460
if (ICOP1(IN) == 0) GO TO 10
IOUT = IOUT+1
BMN(IOUT) = COP(IN)
IBMN(IOUT) = ICOP1(IN)
GO TO 110
10 ICHK = 1
GO TO 110
460 ICHK = 0
NSAVE = ICOP1(IN)
5 continue
do NB=1,NVIRT
  NSIB = MUL(NSM(LN+NB),NSM(NI))
  NSLB = MUL(NSM(LN+NB),LSYM)
  LB = NB-NSYS(NSM(LN+NB))
  INS = NNS(NSIB)
  ILOOP = 0
72 do I=1,INS
    if (INSIN < KBUFF1) GO TO 73
    call dDAFILE(Lu_TiABCI,2,BUFIN,KBUFF1,IAD15)
    INSIN = 0
73  INSIN = INSIN+1
    if (ILOOP == 0) BIAC(I) = BUFIN(INSIN)
    if (ILOOP == 1) BICA(I) = BUFIN(INSIN)
  end do
  ILOOP = ILOOP+1
  if (ILOOP == 1) GO TO 72
  do IT=1,IOUT
    IND = IBMN(IT)
    !ICP1 = mod(IND/IPOW19,8192)
    ICP1 = ibits(IND,19,13)
    INDA = IRC(1)+ICP1
    if (JSYM(INDA) /= NSLB) GO TO 25
    MA = index(INDA)+LB
    !ICP2 = mod(IND/IPOW6,8192)
    !ITYP = mod(IND,64)
    ICP2 = ibits(IND,6,13)
    ITYP = ibits(IND,0,6)
    if (INS == 0) GO TO 25
    COPL = BMN(IT)*C(MA)
    INDB = IRC(ITYP)+ICP2
    FACS = D1
    ICCB = index(INDB)+1
    if (ITYP == 3) GO TO 26
    TERM = DDOT_(INS,C(ICCB),1,BICA,1)
    call DAXPY_(INS,COPL*FACS,BICA,1,S(ICCB),1)
    GO TO 27
26  TERM = DDOT_(INS,C(ICCB),1,BIAC,1)
    call DAXPY_(INS,COPL*FACS,BIAC,1,S(ICCB),1)
27  S(MA) = S(MA)+BMN(IT)*FACS*TERM
25  continue
  end do
end do
if (LEN >= 0) GO TO 100
call DSQ2(C,S,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

return

end subroutine ABCI
