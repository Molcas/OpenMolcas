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

subroutine calcrh1(wrk,wrksize,mapd1,mapd2)
! this routine calcs V1 = V1-V2
!
! mapd1 - direct map of vector 1 (I)
! mapd2 - direct map of vector 2 (I)
!
! N.B. it is assumed, that V1 and V2 are of the same type

#include "wrk.fh"
integer mapd1(0:512,1:6)
integer mapd2(0:512,1:6)
! help variables
integer poss1, poss2, length, ii

!1 calc length
ii = mapd1(0,5)
length = mapd1(ii,1)+mapd1(ii,2)-mapd1(1,1)

!2 realize substract
if (length > 0) then
  poss1 = mapd1(1,1)
  poss2 = mapd2(1,1)
  do ii=0,length-1
    wrk(poss1+ii) = wrk(poss1+ii)-wrk(poss2+ii)
  end do
end if

return

end subroutine calcrh1
