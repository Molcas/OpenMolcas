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

subroutine iRBuf(Array,nArray,Copy)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp), intent(in) :: nArray
integer(kind=iwp), intent(inout) :: Array(nArray)
logical(kind=iwp), intent(in) :: Copy

call idRBuf(Array,nArray/RtoI,Copy)

return

! This is to allow type punning without an explicit interface
contains

subroutine idRBuf(Array,nArray,Copy)

  integer(kind=iwp), intent(in) :: nArray
  logical(kind=iwp), intent(in) :: Copy
  integer(kind=iwp), target, intent(inout) :: Array(nArray)
  real(kind=wp), pointer :: dArray(:)

  call c_f_pointer(c_loc(Array(1)),dArray,[nArray])
  call dRBuf(dArray,nArray,Copy)
  nullify(dArray)

end subroutine idRBuf

end subroutine iRBuf
