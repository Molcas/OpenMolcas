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

subroutine dRdMCK(rc,Option,InLab,iComp,dData,iSymLab)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp), intent(in) :: Option, iComp, iSymLab
character(len=*), intent(inout) :: InLab
real(kind=wp), intent(_OUT_) :: dData(*)

call dRdMCK_Internal(dData)

! This is to allow type punning without an explicit interface
contains

subroutine dRdMCK_Internal(dData)

  real(kind=wp), target, intent(_OUT_) :: dData(*)
  integer(kind=iwp), pointer :: iData(:)

  call c_f_pointer(c_loc(dData),iData,[1])
  call RdMCK(rc,Option,InLab,iComp,iData,iSymLab)
  nullify(iData)

  return

end subroutine dRdMCK_Internal

end subroutine dRdMCK
