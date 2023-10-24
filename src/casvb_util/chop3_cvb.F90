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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine chop3_cvb()

use casvb_global, only: aikcof, bikcof, cikcof, ifnss1, ifnss2, ikcoff, kbasis, kbasiscvb, nalf, nbet, ndetvbs, nel, release
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iretval1, iretval2, kmost, mxdetvb, mxfns, need

if (release(3)) then
  call mma_deallocate(aikcof)
  nullify(bikcof)
  if (allocated(cikcof)) call mma_deallocate(cikcof)
  call mma_deallocate(ikcoff)
  call mma_deallocate(ifnss1)
  call mma_deallocate(ifnss2)
  call mma_deallocate(ndetvbs)
end if
release(3) = .true.
release(4) = .false.

call icomb_cvb(nel,nbet,iretval1)
call icomb_cvb(nel,nbet-1,iretval2)
mxfns = iretval1-iretval2
if (kbasis == 5) call icomb_cvb(nel,nalf,mxfns)
call icomb_cvb(nel,nalf,mxdetvb)
if (((kbasis > 2) .and. (kbasis /= 6)) .or. ((kbasiscvb > 2) .and. (kbasiscvb /= 6))) then
  kmost = 3
else if ((kbasis <= 2) .or. (kbasiscvb <= 2)) then
  kmost = 1
else
  kmost = 6
end if
call bspset_cvb(kmost,1,need)
if (kmost == 3) then
  call mma_allocate(aikcof,[0,need],label='aikcof')
  call mma_allocate(cikcof,[0,need],label='cikcof')
  bikcof => cikcof
else if (kmost == 1) then
  call mma_allocate(aikcof,[0,need],label='aikcof')
  bikcof => aikcof
else
  call mma_allocate(aikcof,[0,0],label='aikcof')
  bikcof => aikcof
end if
! Flag AIKCOF/BIKCOF as unset:
aikcof(0) = Zero
bikcof(0) = Zero

call mma_allocate(ikcoff,[0,nel],[0,nel],[0,nel],label='ikcoff')
call mma_allocate(ifnss1,[0,nel],[0,nel],label='ifnss1')
call mma_allocate(ifnss2,[0,nel],[0,nel],label='ifnss2')
call mma_allocate(ndetvbs,[0,nel],[0,nel],label='ndetvbs')
call bspset_cvb(kbasiscvb,2,need)

return

end subroutine chop3_cvb
