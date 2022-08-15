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

subroutine DMSMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: MemDer, MmEFP, nOrder

Mem = 0
nHer = 0
call EFMmP(nOrder,MmEFP,la,lb+1,lr-1)
Mem = max(Mem,MmEFP)
nHer = max(nHer,nOrder)
call EFMmP(nOrder,MmEFP,la,lb,lr-1)
Mem = max(Mem,MmEFP)
nHer = max(nHer,nOrder)

! Add a scratch area for intermediate integrals

MemDer = 3*(nTri_Elem1(la)*nTri_Elem1(lb+1)+nTri_Elem1(la)*nTri_Elem1(lb))
Mem = Mem+MemDer
Mem = Mem+nTri_Elem1(la)*nTri_Elem1(lb)*9

return

end subroutine DMSMem
