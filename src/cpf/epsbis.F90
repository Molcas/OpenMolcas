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

subroutine EPSBIS(JSY,INDX,C,W,EPB)

use cpf_global, only: ICPF, INCPF, IPRINT, IRC, ISDCI, LSYM, NNS, NVIR
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*)
real(kind=wp), intent(in) :: C(*), W(*)
real(kind=wp), intent(_OUT_) :: EPB(*)
integer(kind=iwp) :: I, IIN, INUM, IP, IST, NS1, NSIL
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

EPB(1:IRC(4)) = Zero
if ((ICPF == 1) .or. (ISDCI == 1) .or. (INCPF == 1)) return

! VALENCE

IP = IRC(1)
do I=1,IP
  EPB(I) = C(I)*W(I)
end do

! SINGLES

IP = IRC(2)-IRC(1)
IIN = IRC(1)
do I=1,IP
  !FUE IND = IND+1
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  EPB(IIN+I) = DDOT_(INUM,C(IST),1,W(IST),1)
end do

! DOUBLES

IP = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IP
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  EPB(IIN+I) = DDOT_(INUM,C(IST),1,W(IST),1)
end do

IP = IRC(4)
if (IPRINT > 5) write(u6,998) (EPB(I),I=1,IP)

return

998 format(6X,'EPB ',5F10.6)

end subroutine EPSBIS
