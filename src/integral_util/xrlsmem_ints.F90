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

subroutine xRlsMem_Ints()

use k2_arrays, only: XMem

implicit none

if (Xmem) then
  !write(u6,*)
  !write(u6,*) 'xRlsMem_Ints: External allocation deactivate!'
  !write(u6,*)
  XMem = .false.
  call RlsMem_Ints()
  !write(u6,*) 'xRlsMem_Ints: External allocation released!'
!else
!  write(u6,*)
!  write(u6,*) 'xRlsMem_Ints: No external scratch handling to deactivate!'
!  write(u6,*)
end if

end subroutine xRlsMem_Ints
