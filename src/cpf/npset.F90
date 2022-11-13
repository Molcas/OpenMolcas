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

subroutine NPSET(JSY,INDX,C,TPQ,ENP,T,S,W,EPP,ICASE)

use cpf_global, only: IADDP, ICPF, IDENS, INCPF, IPRINT, IRC, ISDCI, ITPUL, JSC, LSYM, Lu_CI, NCONF, NNS, NVIR
use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*), ICASE(*)
real(kind=wp), intent(inout) :: C(*), ENP(*)
real(kind=wp), intent(_OUT_) :: TPQ(*), T(*), S(*), W(*), EPP(*)
integer(kind=iwp) :: I, IAD, IIN, IND, INUM, IP, IQ, IST, NS1, NSIL
real(kind=wp) :: EMPI
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

if (IDENS /= 1) then

  if (ITPUL == 1) then
    IAD = 0
    IADDP(1) = 0
    call dDAFILE(Lu_CI,1,C,NCONF,IAD)
    IADDP(2) = IAD
  end if

  ! VALENCE

  IQ = IRC(1)
  T(1:IQ) = C(1:IQ)**2

  ! SINGLES

  IQ = IRC(2)-IRC(1)
  IND = IRC(1)
  IIN = IRC(1)
  do I=1,IQ
    IND = IND+1
    NS1 = JSUNP(JSY,IIN+I)
    NSIL = MUL(NS1,LSYM)
    INUM = NVIR(NSIL)
    IST = INDX(IIN+I)+1
    T(IND) = DDOT_(INUM,C(IST),1,C(IST),1)
  end do

  ! DOUBLES

  IQ = IRC(4)-IRC(2)
  IIN = IRC(2)
  do I=1,IQ
    IND = IND+1
    NS1 = JSUNP(JSY,IIN+I)
    NSIL = MUL(NS1,LSYM)
    INUM = NNS(NSIL)
    IST = INDX(IIN+I)+1
    T(IND) = DDOT_(INUM,C(IST),1,C(IST),1)
  end do
  IP = IRC(4)
  do I=1,IP
    call TPQSET(ICASE,TPQ,I)
    ENP(I) = DDOT_(IP,TPQ,1,T,1)+One
  end do
  IP = IRC(4)
  if (IPRINT > 5) write(u6,12) (ENP(I),I=1,IP)

end if

! VALENCE

IQ = IRC(1)
do I=1,IQ
  if (IDENS == 0) then
    EMPI = One/sqrt(ENP(I))
  else
    EMPI = sqrt(ENP(I))
  end if
  C(I) = C(I)*EMPI
end do

! SINGLES

IQ = IRC(2)-IRC(1)
IIN = IRC(1)
do I=1,IQ
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  if (IDENS == 0) then
    EMPI = One/sqrt(ENP(IIN+I))
  else
    EMPI = sqrt(ENP(IIN+I))
  end if
  C(IST:IST+INUM-1) = EMPI*C(IST:IST+INUM-1)
end do

! DOUBLES

IQ = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IQ
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  if (IDENS == 0) then
    EMPI = One/sqrt(ENP(IIN+I))
  else
    EMPI = sqrt(ENP(IIN+I))
  end if
  C(IST:IST+INUM-1) = EMPI*C(IST:IST+INUM-1)
end do
if (IPRINT >= 15) write(u6,13) (C(I),I=1,NCONF)

if (IDENS /= 1) then
  EPP(1:IRC(4)) = Zero
  S(1:JSC(4)) = Zero
  if ((ICPF /= 1) .and. (ISDCI /= 1) .and. (INCPF /= 1)) W(1:JSC(4)) = Zero
end if

return

12 format(6X,'ENP ',5F14.8)
13 format(6X,'C(NP)',5F10.6)

end subroutine NPSET
