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
function len_trim_cvb(a)
! Length of string excluding trailing blanks

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: len_trim_cvb
character(len=*) :: a
integer(kind=iwp) :: i

do i=len(a),1,-1
  if (a(i:i) /= ' ') exit
end do
len_trim_cvb = i

return

end function len_trim_cvb
