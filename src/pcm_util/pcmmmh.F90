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

subroutine PCMMmh( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: iAng(4)

iAng(1) = la
iAng(2) = lb
iAng(3) = lr
iAng(4) = 0
call MemRg2(iAng,nHer,Mem,2)
Mem = Mem+2+nTri_Elem1(la)*nTri_Elem1(lb)*nTri_Elem1(lr)

return

end subroutine PCMMmh
