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

!pgi$g opt=1
subroutine MABCI(JSY,INDX,C,S,BMN,IBMN,BIAC,BICA,BUFIN,W,THET,ENP,NII)

use Constants, only: One, Two, Half
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp) :: JSY(*), INDX(*), IBMN(*), NII
real(kind=wp) :: C(*), S(*), BMN(*), BIAC(*), BICA(*), BUFIN(*), W(*), THET(NII,NII), ENP(*)
#include "cpfmcpf.fh"
#include "files_cpf.fh"
integer(kind=iwp) :: I, IAD15, ICCB, ICHK, ICP1, ICP2, IIN, ILEN, ILOOP, IND, INDA, INDB, INS, INSIN, INUM, IOUT, IT, ITYP, LB, &
                     MA, NB, NI, NSAVE, NSIB, NSLB
real(kind=wp) :: COPL, ENPQ, FACS, FACW, FACWA, FACWB, TERM, XXX
integer(kind=iwp), external :: JSUNP_CPF
real(kind=r8), external :: DDOT_
!parameter(IPOW6=2**6,IPOW13=2**13,IPOW19=2**19)
! Statement function
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP_CPF(JSY,L)

INUM = IRC(4)-IRC(3)
call MPSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
ICHK = 0
INSIN = KBUFF1
IAD15 = IADABCI
IADD10 = IAD10(4)
call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
IIN = 2
NSAVE = ICOP1(IIN)
100 NI = NSAVE
IOUT = 0
110 IIN = IIN+1
if (IIN <= ILEN) GO TO 15
call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
if (ILEN <= 0) GO TO 5
IIN = 1
15 if (ICHK /= 0) GO TO 460
if (ICOP1(IIN) == 0) GO TO 10
IOUT = IOUT+1
BMN(IOUT) = COP(IIN)
IBMN(IOUT) = ICOP1(IIN)
GO TO 110
10 ICHK = 1
GO TO 110
460 ICHK = 0
NSAVE = ICOP1(IIN)
5 continue
do NB=1,NVIRT
  NSIB = MUL(NSM(LN+NB),NSM(NI))
  NSLB = MUL(NSM(LN+NB),LSYM)
  LB = NB-NSYS(NSM(LN+NB))
  INS = NNS(NSIB)
  ILOOP = 0
72 continue
  do I=1,INS
    if (INSIN >= KBUFF1) then
      call dDAFILE(Lu_TiABCI,2,BUFIN,KBUFF1,IAD15)
      INSIN = 0
    end if
    INSIN = INSIN+1
    if (ILOOP == 0) BIAC(I) = BUFIN(INSIN)
    if (ILOOP == 1) BICA(I) = BUFIN(INSIN)
  end do
  ILOOP = ILOOP+1
  if (ILOOP == 1) GO TO 72
  do IT=1,IOUT
    IND = IBMN(IT)
    !PAM97 ICP1 = iand(ishft(IND,-19),8191)
    !ICP1 = mod(IND/IPOW19,IPOW13)
    ICP1 = ibits(IND,19,13)
    INDA = IRC(1)+ICP1
    if (JSYM(INDA) /= NSLB) GO TO 25
    MA = INDX(INDA)+LB
    !PAM97 ICP2 = iand(ishft(IND,-6),8191)
    !PAM97 ITYP = iand(IND,63)
    !ICP2 = mod(IND/IPOW6,IPOW13)
    !ITYP = mod(IND,IPOW6)
    ICP2 = ibits(IND,6,13)
    ITYP = ibits(IND,0,6)
    if (INS == 0) GO TO 25
    COPL = BMN(IT)*C(MA)
    INDB = IRC(ITYP)+ICP2
    XXX = THET(INDA,INDB)*Half
    ENPQ = (One-XXX)*(ENP(INDA)+ENP(INDB)-One)+XXX
    FACS = sqrt(ENP(INDA))*sqrt(ENP(INDB))/ENPQ
    FACW = FACS*(Two-THET(INDA,INDB))/ENPQ
    FACWA = FACW*ENP(INDA)-FACS
    FACWB = FACW*ENP(INDB)-FACS
    ICCB = INDX(INDB)+1
    if (ITYP == 3) then
      TERM = DDOT_(INS,C(ICCB),1,BIAC,1)
      call DAXPY_(INS,COPL*FACS,BIAC,1,S(ICCB),1)
      call DAXPY_(INS,COPL*FACWB,BIAC,1,W(ICCB),1)
    else
      TERM = DDOT_(INS,C(ICCB),1,BICA,1)
      call DAXPY_(INS,COPL*FACS,BICA,1,S(ICCB),1)
      call DAXPY_(INS,COPL*FACWB,BICA,1,W(ICCB),1)
    end if
    S(MA) = S(MA)+BMN(IT)*FACS*TERM
    W(MA) = W(MA)+BMN(IT)*FACWA*TERM
25  continue
  end do
end do
if (ILEN >= 0) GO TO 100
call MDSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

return

end subroutine MABCI
