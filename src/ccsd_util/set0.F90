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

subroutine set0(wrk,wrksize,mapd,mapi)
! this routine vanishes given mediate

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapd(0:512,6), mapi(8,8,8)
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: ii, length, poss0

!1 def poss0, legth
poss0 = mapd(1,1)
ii = mapd(0,5)
length = mapd(ii,1)+mapd(ii,2)-poss0

!2 set appropriate wrk = 0
call mv0zero(length,length,wrk(poss0))

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(mapi)

end subroutine set0
