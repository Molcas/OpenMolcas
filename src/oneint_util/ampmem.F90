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
! Copyright (C) 1996, Per Ake Malmqvist                                *
!               1996, Roland Lindh                                     *
!***********************************************************************

subroutine AMPMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: Mem1, Mem2, Mem3, nOrder

#include "macros.fh"
unused_var(lr)

! Mem1: Workspace for MltPrm.
! Mem2: Tables Tpp,Tp,T0,Tm, and Tmm.
! Mem3: Result from AMPr.
Mem1 = 0
call MltMmP(nOrder,Mem,la,lb+2,2)
Mem1 = max(Mem1,Mem)
Mem2 = 6*nTri_Elem1(la)*nTri_Elem1(lb+2)
nHer = nOrder
call MltMmP(nOrder,Mem,la,lb+1,1)
Mem1 = max(Mem1,Mem)
Mem2 = Mem2+3*nTri_Elem1(la)*nTri_Elem1(lb+1)
call MltMmP(nOrder,Mem,la,lb,2)
Mem1 = max(Mem1,Mem)
Mem2 = Mem2+6*nTri_Elem1(la)*nTri_Elem1(lb)
if (lb >= 1) then
  call MltMmP(nOrder,Mem,la,lb-1,1)
  Mem1 = max(Mem1,Mem)
  Mem2 = Mem2+3*nTri_Elem1(la)*nTri_Elem1(lb-1)
  if (lb >= 2) then
    call MltMmP(nOrder,Mem,la,lb-2,2)
    Mem1 = max(Mem1,Mem)
    Mem2 = Mem2+6*nTri_Elem1(la)*nTri_Elem1(lb-2)
  end if
end if
Mem3 = 6*nTri_Elem1(la)*nTri_Elem1(lb)

Mem = Mem1+Mem2+Mem3+1

return

end subroutine AMPMem
