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

subroutine pXpMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: MemPX, nOrder

Mem = 0
nHer = 0
call PXMem(nOrder,MemPX,la,lb+1,lr-1)
Mem = max(Mem,MemPX)
nHer = max(nHer,nOrder)
if (lb > 0) then
  call PXMem(nOrder,MemPX,la,lb-1,lr-1)
  Mem = max(Mem,MemPX)
  nHer = max(nHer,nOrder)
end if

return

end subroutine pXpMem
