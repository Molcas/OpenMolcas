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

use cpf_global, only: IADDP, ICPF, IDENS, INCPF, IPRINT, IRC, ISDCI, ITPUL, JSC, LSYM, Lu_CI, MUL, NCONF, NNS, NVIR
use Constants, only: One
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp) :: JSY(*), INDX(*), ICASE(*)
real(kind=wp) :: C(*), TPQ(*), ENP(*), T(*), S(*), W(*), EPP(*)
integer(kind=iwp) :: I, IAD, IIN, IND, INUM, IP, IQ, IST, NS1, NSIL
real(kind=wp) :: EMPI
integer(kind=iwp), external :: JSUNP_CPF
real(kind=r8), external :: DDOT_

if (IDENS /= 1) then

  if (ITPUL == 1) then
    IAD = 0
    IADDP(1) = 0
    call dDAFILE(Lu_CI,1,C,NCONF,IAD)
    IADDP(2) = IAD
  end if

  ! VALENCE

  IQ = IRC(1)
  do I=1,IQ
    T(I) = C(I)*C(I)
  end do

  ! SINGLES

  IQ = IRC(2)-IRC(1)
  IND = IRC(1)
  IIN = IRC(1)
  do I=1,IQ
    IND = IND+1
    NS1 = JSUNP_CPF(JSY,IIN+I)
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
    NS1 = JSUNP_CPF(JSY,IIN+I)
    NSIL = MUL(NS1,LSYM)
    INUM = NNS(NSIL)
    IST = INDX(IIN+I)+1
    T(IND) = DDOT_(INUM,C(IST),1,C(IST),1)
  end do
  IP = IRC(4)
  do I=1,IP
    call TPQSET(ICASE,TPQ,I)
    ENP(I) = DDOT_(IP,TPQ,1,T,1)
    ENP(I) = ENP(I)+One
  end do
  IP = IRC(4)
  if (IPRINT > 5) write(u6,12) (ENP(I),I=1,IP)

end if

! VALENCE

IQ = IRC(1)
do I=1,IQ
  if (IDENS == 0) EMPI = One/sqrt(ENP(I))
  if (IDENS == 1) EMPI = sqrt(ENP(I))
  C(I) = C(I)*EMPI
end do

! SINGLES

IQ = IRC(2)-IRC(1)
IIN = IRC(1)
do I=1,IQ
  NS1 = JSUNP_CPF(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  if (IDENS == 0) EMPI = One/sqrt(ENP(IIN+I))
  if (IDENS == 1) EMPI = sqrt(ENP(IIN+I))
  call VSMUL(C(IST),1,EMPI,C(IST),1,INUM)
end do

! DOUBLES

IQ = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IQ
  NS1 = JSUNP_CPF(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  if (IDENS == 0) EMPI = One/sqrt(ENP(IIN+I))
  if (IDENS == 1) EMPI = sqrt(ENP(IIN+I))
  call VSMUL(C(IST),1,EMPI,C(IST),1,INUM)
end do
if (IPRINT >= 15) write(u6,13) (C(I),I=1,NCONF)
if (IDENS == 1) return

call SETZ(EPP,IRC(4))
call SETZ(S,JSC(4))
if ((ICPF /= 1) .and. (ISDCI /= 1) .and. (INCPF /= 1)) call SETZ(W,JSC(4))

return

12 format(6X,'ENP ',5F14.8)
13 format(6X,'C(NP)',5F10.6)

end subroutine NPSET
