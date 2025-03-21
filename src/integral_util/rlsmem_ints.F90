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

subroutine rlsmem_ints()

use k2_arrays, only: Sew_Scr, XMem
use stdalloc, only: mma_deallocate

implicit none

if (.not. XMem) then
  call mma_deallocate(Sew_Scr,safe='*')
!else
!  write(u6,*) 'RlsMem_Ints: External scratch handling active!'
end if

return

end subroutine rlsmem_ints
