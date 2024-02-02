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

implicit none
integer n, a(n), temp, k, l

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
