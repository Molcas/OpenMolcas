!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine PVMem(nHer,Mem,la,lb,lr,KrnMem)

use Integral_interfaces, only: int_mem
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: nHer, Mem
integer(kind=iwp), intent(in) :: la, lb, lr
procedure(int_mem) :: KrnMem
integer(kind=iwp) :: MemNA1, MemNA2

call KrnMem(nHer,MemNA1,la+1,lb,lr-1)

if (la /= 0) then
  call KrnMem(nHer,MemNA2,la-1,lb,lr-1)
else
  MemNA2 = 0
end if

Mem = max(MemNA1,MemNA2)

return

end subroutine PVMem
