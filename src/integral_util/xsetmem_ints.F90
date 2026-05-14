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

subroutine xsetmem_ints(mem)

use k2_arrays, only: Sew_Scr, XMem
use stdalloc, only: mma_allocate, mma_maxDBLE
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: mem
integer(kind=iwp) :: mem_, MemMax

if (XMem) then
  call WarningMessage(2,'External handling of scratch already active!')
  call Abend()
end if
Mem_ = Mem
! Avoid using up all available memory
call mma_maxDBLE(MemMax)
if ((MemMax-Mem_ < 8000) .and. (Mem_ > 8000)) Mem_ = Mem_-8000
!write(u6,*) 'xsetmem_ints: External allocate:',Mem_
call mma_allocate(Sew_Scr,Mem_,Label='Sew_Scr')
XMem = .true.
!call mma_MaxDBLE(nu)
!write(u6,*) 'xsetmem_ints: External allocate left to allocate:',nu

end subroutine xsetmem_ints
