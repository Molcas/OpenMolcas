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

subroutine MABCI(JSY,INDX,C,S,BMN,IBMN,BIAC,BICA,BUFIN,W,THET,ENP,NII)

use cpf_global, only: IADABCI, IRC, KBUFF1, LN, LSYM, Lu_CIGuga, Lu_TiABCI, NDIAG, NNS, NSM, NSYS, NVIRT, SQ2
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*), NII
real(kind=wp), intent(_OUT_) :: C(*), S(*), W(*)
real(kind=wp), intent(inout) :: BMN(*), BIAC(*), BICA(*), BUFIN(*)
integer(kind=iwp), intent(_OUT_) :: IBMN(*)
real(kind=wp), intent(in) :: THET(NII,NII), ENP(*)
integer(kind=iwp) :: I, IAD15, IADD10, ICCB, ICHK, ICP1, ICP2, IIN, ILEN, ILOOP, IND, INDA, INDB, INS, INSIN, INUM, IOUT, IT, &
                     ITYP, LB, MA, NB, NI, NSAVE, NSIB, NSLB
real(kind=wp) :: COPL, ENPQ, FACS, FACW, FACWA, FACWB, TERM, XXX
logical(kind=iwp) :: Skip
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

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
do
  NI = NSAVE
  IOUT = 0
  Skip = .false.
  do
    IIN = IIN+1
    if (IIN > ILEN) then
      call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      ILEN = ICOP1(nCOP+1)
      if (ILEN <= 0) then
        Skip = .true.
        exit
      end if
      IIN = 1
    end if
    if (ICHK /= 0) exit
    if (ICOP1(IIN) /= 0) then
      IOUT = IOUT+1
      BMN(IOUT) = COP(IIN)
      IBMN(IOUT) = ICOP1(IIN)
    else
      ICHK = 1
    end if
  end do
  if (.not. Skip) then
    ICHK = 0
    NSAVE = ICOP1(IIN)
  end if
  do NB=1,NVIRT
    NSIB = MUL(NSM(LN+NB),NSM(NI))
    NSLB = MUL(NSM(LN+NB),LSYM)
    LB = NB-NSYS(NSM(LN+NB))
    INS = NNS(NSIB)
    ILOOP = 0
    do
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
      if (ILOOP /= 1) exit
    end do
    do IT=1,IOUT
      IND = IBMN(IT)
      ICP1 = ibits(IND,19,13)
      INDA = IRC(1)+ICP1
      if (JSUNP(JSY,INDA) /= NSLB) cycle
      MA = INDX(INDA)+LB
      ICP2 = ibits(IND,6,13)
      ITYP = ibits(IND,0,6)
      if (INS == 0) cycle
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
        S(ICCB:ICCB+INS-1) = S(ICCB:ICCB+INS-1)+COPL*FACS*BIAC(1:INS)
        W(ICCB:ICCB+INS-1) = W(ICCB:ICCB+INS-1)+COPL*FACWB*BIAC(1:INS)
      else
        TERM = DDOT_(INS,C(ICCB),1,BICA,1)
        S(ICCB:ICCB+INS-1) = S(ICCB:ICCB+INS-1)+COPL*FACS*BICA(1:INS)
        W(ICCB:ICCB+INS-1) = W(ICCB:ICCB+INS-1)+COPL*FACWB*BICA(1:INS)
      end if
      S(MA) = S(MA)+BMN(IT)*FACS*TERM
      W(MA) = W(MA)+BMN(IT)*FACWA*TERM
    end do
  end do
  if (ILEN < 0) exit
end do
call MDSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

return

end subroutine MABCI
