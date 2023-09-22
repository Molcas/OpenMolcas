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

!IFG trivial
subroutine imove_cvb(iv1,iv2,n)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: n, iv1(n), iv2(n)
integer(kind=iwp) :: i

do i=1,n
  iv2(i) = iv1(i)
end do

return

end subroutine imove_cvb
