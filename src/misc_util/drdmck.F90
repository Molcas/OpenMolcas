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

implicit integer(A-Z)
character*(*) InLab
real*8 dData(*)

call dRdMCK_Internal(dData)

! This is to allow type punning without an explicit interface
contains

subroutine dRdMCK_Internal(dData)

  use iso_c_binding

  real*8, target :: dData(*)
  integer, pointer :: iData(:)

  call c_f_pointer(c_loc(dData),iData,[1])
  call RdMCK(rc,Option,InLab,iComp,iData,iSymLab)
  nullify(iData)

  return

end subroutine dRdMCK_Internal

end subroutine dRdMCK
