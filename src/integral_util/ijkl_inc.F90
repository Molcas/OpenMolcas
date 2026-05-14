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

subroutine ijkl_inc(i,j,k,l)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: i, j, k, l

l = l+1
if (((i == k) .and. (l > j)) .or. ((i /= k) .and. (l > k))) then
  l = 1
  k = k+1
  if (k > i) then
    k = 1
    j = j+1
    if (j > i) then
      j = 1
      i = i+1
    end if
  end if
end if

return

end subroutine ijkl_inc
