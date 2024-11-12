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

subroutine pa_sort(a,n)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(inout) :: a(n)
integer(kind=iwp) :: k, l, temp

do k=1,n-1
  do l=k+1,n
    if (a(k) > a(l)) then
      temp = a(k)
      a(k) = a(l)
      a(l) = temp
    end if
  end do
end do

return

end subroutine pa_sort
