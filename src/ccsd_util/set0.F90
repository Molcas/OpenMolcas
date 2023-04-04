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

subroutine set0(wrk,wrksize,map)
! this routine vanishes given mediate

use ccsd_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize
real(kind=wp) :: wrk(wrksize)
type(Map_Type) :: map
integer(kind=iwp) :: ii, length, poss0

!1 def poss0, legth
poss0 = map%d(1,1)
ii = map%d(0,5)
length = map%d(ii,1)+map%d(ii,2)-poss0

!2 set appropriate wrk = 0
call mv0zero(length,length,wrk(poss0))

return

end subroutine set0
