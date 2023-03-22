!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine GEN_INT_DCCD(rc,ipq1,Xint)

use Index_Functions, only: nTri_Elem
use GetInt_mod, only: nRS, Vec2, NumV
use GetInt_mod, only: lists, I, hash_table
use TwoDat, only: rcTwo
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: ipq1
real(kind=wp), intent(_OUT_) :: Xint(*)
integer(kind=iwp) :: iR, iR_, iRS, iS, iS_, J
real(kind=wp) :: Temp

Xint(1:nRS) = Zero

do iR_=lists(3,I),lists(4,I)
  iR = hash_table(iR_)
  do iS_=lists(3,I),iR_
    iS = hash_table(iS_)
    iRS = nTri_Elem(iR-1)+iS
    Temp = Zero
    do J=1,NumV
      Temp = Temp+Vec2(iRS,J)*Vec2(ipq1,J)
    end do
    XInt(iRS) = XInt(iRS)+Temp
  end do
end do

rc = rcTwo%good

return

end subroutine GEN_INT_DCCD
