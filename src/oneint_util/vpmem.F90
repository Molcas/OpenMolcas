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

subroutine VPMem( &
#                define _CALLING_
#                include "mem_interface.fh"
                )

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: iAngV(4), MemNA1, MemNA2, nFlop, nMem

call mHrr(la,lb+1,nFlop,nMem)

nHer = (la+lb+1+lr+2)/2
iAngV(1) = la
iAngV(2) = lb+1
iAngV(3) = 0
iAngV(4) = 0
call MemRys(iAngV,MemNA1)

MemNA1 = max(nMem,MemNA1)

if (lb /= 0) then
  call mHrr(la,lb-1,nFlop,nMem)

  nHer = (la+lb-1+lr+2)/2
  iAngV(1) = la
  iAngV(2) = lb-1
  iAngV(3) = 0
  iAngV(4) = 0
  call MemRys(iAngV,MemNA2)

  MemNA2 = max(nMem,MemNA2)
else
  MemNA2 = 0
end if

Mem = max(MemNA1,MemNA2)

Mem = Mem+1

Mem = Mem+nTri_Elem1(la)*nTri_Elem1(lb+1)
if (lb /= 0) Mem = Mem+nTri_Elem1(la)*nTri_Elem1(lb-1)

return

end subroutine VPMem
