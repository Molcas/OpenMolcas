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

subroutine EPSPRIM(JSY,INDX,C,S,EPP)

use cpf_global, only: IPRINT, IRC, LSYM, NNS, NVIR
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*)
real(kind=wp), intent(in) :: C(*), S(*)
real(kind=wp), intent(inout) :: EPP(*)
integer(kind=iwp) :: I, IIN, INUM, IP, IST, NS1, NSIL
real(kind=wp) :: T
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

! VALENCE

IP = IRC(1)
do I=1,IP
  EPP(I) = EPP(I)+C(I)*S(I)
end do

! SINGLES

IP = IRC(2)-IRC(1)
IIN = IRC(1)
do I=1,IP
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  T = DDOT_(INUM,C(IST),1,S(IST),1)
  EPP(IIN+I) = EPP(IIN+I)+T
end do

! DOUBLES

IP = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IP
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  T = DDOT_(INUM,C(IST),1,S(IST),1)
  EPP(IIN+I) = EPP(IIN+I)+T
end do

IP = IRC(4)
if (IPRINT > 5) write(u6,998) (EPP(I),I=1,IP)

return

998 format(6X,'EPP ',5F10.6)

end subroutine EPSPRIM
