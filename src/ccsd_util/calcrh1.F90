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

subroutine calcrh1(wrk,wrksize,v1,v2)
! this routine calcs V1 = V1-V2
!
! v1 - map type of vector 1 (I)
! v2 - map type of vector 2 (I)
!
! N.B. it is assumed, that V1 and V2 are of the same type

use ccsd_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize
real(kind=wp) :: wrk(wrksize)
type(Map_Type) :: v1, v2
integer(kind=iwp) :: ii, length, poss1, poss2

!1 calc length
ii = v1%d(0,5)
length = v1%d(ii,1)+v1%d(ii,2)-v1%d(1,1)

!2 realize substract
if (length > 0) then
  poss1 = v1%d(1,1)
  poss2 = v2%d(1,1)
  do ii=0,length-1
    wrk(poss1+ii) = wrk(poss1+ii)-wrk(poss2+ii)
  end do
end if

return

end subroutine calcrh1
