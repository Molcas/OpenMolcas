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

subroutine WrOne(rc,Option,InLab,Comp,rData,SymLab)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: rc, Option, Comp, SymLab
character(len=*) :: InLab
real(kind=wp) :: rData(*)

call WrOne_Internal(rData)

! This is to allow type punning without an explicit interface
contains

subroutine WrOne_Internal(rData)

  real(kind=wp), target :: rData(*)
  integer(kind=iwp), pointer :: iData(:)

  call c_f_pointer(c_loc(rData(1)),iData,[1])
  call iWrOne(rc,Option,InLab,Comp,iData,SymLab)
  nullify(iData)

  return

end subroutine WrOne_Internal

end subroutine WrOne
