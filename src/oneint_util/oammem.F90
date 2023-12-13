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

subroutine OAMMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: MmMltP, nOrder

call MltMmP(nOrder,Mem,la,lb+1,lr-1)
nHer = nOrder
if (lb > 0) then
  call MltMmP(nOrder,MmMltP,la,lb-1,lr-1)
  Mem = max(Mem,MmMltP)+nTri_Elem1(la)*nTri_Elem1(lb-1)*3
end if
Mem = Mem+1+nTri_Elem1(la)*nTri_Elem1(lb+1)*3
Mem = Mem+nTri_Elem1(la)*nTri_Elem1(lb)*3

return

end subroutine OAMMem
