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

subroutine cWrMCK(rc,Option,InLab,iComp,cData,iSymLab)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp), intent(in) :: Option, iComp, iSymLab
character(len=*), intent(in) :: InLab
character, intent(in) :: cData(*)

call cWrMCK_Internal(cData)

! This is to allow type punning without an explicit interface
contains

subroutine cWrMCK_Internal(cData)

  character, target, intent(in) :: cData(*)
  integer(kind=iwp), pointer :: iData(:)

  call c_f_pointer(c_loc(cData(1)),iData,[1])
  call WrMCK(rc,Option,InLab,iComp,iData,iSymLab)
  nullify(iData)

  return

end subroutine cWrMCK_Internal

end subroutine cWrMCK
