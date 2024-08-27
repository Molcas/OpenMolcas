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

function iGet(A,n)

use iso_c_binding

implicit none
integer iGet
integer n
real*8, target :: A(*)
integer, pointer :: iA(:)

call c_f_pointer(c_loc(A),iA,[n])
iGet = iA(n)
nullify(iA)

return

end function iGet
