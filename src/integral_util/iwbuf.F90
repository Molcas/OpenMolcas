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

subroutine iWBuf(Array,nArray)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp), intent(in) :: nArray, Array(nArray)

call idWBuf(Array,nArray/RtoI)

return

! This is to allow type punning without an explicit interface
contains

subroutine idWBuf(Array,nArray)

  integer(kind=iwp), intent(in) :: nArray
  integer(kind=iwp), target, intent(in) :: Array(nArray)
  real(kind=wp), pointer :: dArray(:)

  call c_f_pointer(c_loc(Array(1)),dArray,[nArray])
  call dWBuf(dArray,nArray)
  nullify(dArray)

end subroutine idWBuf

end subroutine iWBuf
