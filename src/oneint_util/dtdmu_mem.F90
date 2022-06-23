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

subroutine dTdmu_mem( &
#                    define _CALLING_
#                    include "mem_interface.fh"
                    )

use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: MemDer, MmEFP, nOrder
! Statement function for Cartesian index
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

Mem = 0
nHer = 0
call EFMmP(nOrder,MmEFP,la,lb+1,lr)
Mem = max(Mem,MmEFP)
nHer = max(nHer,nOrder)
if (lb >= 1) then
  call EFMmP(nOrder,MmEFP,la,lb-1,lr)
  Mem = max(Mem,MmEFP)
  nHer = max(nHer,nOrder)
end if

! Add a scratch area for intermediate integrals

MemDer = 3*nElem(la)*nElem(lb+1)
if (lb >= 1) MemDer = MemDer+3*nElem(la)*nElem(lb-1)
Mem = Mem+MemDer+1
Mem = Mem+nElem(la)*nElem(lb)*3

return

end subroutine dTdmu_mem
