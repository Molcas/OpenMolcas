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

function C2R8(CBuf)

real*8 C2R8
character CBuf(*)

C2R8 = C2R8_Internal(CBuf)

! This is to allow type punning without an explicit interface
contains

real*8 function C2R8_Internal(CBuf)

  use iso_c_binding

  character, target :: CBuf(*)
  real*8, pointer :: Buf

  call c_f_pointer(c_loc(CBuf(1)),Buf)
  C2R8_Internal = Buf
  nullify(Buf)

  return

end function C2R8_Internal

end function C2R8
