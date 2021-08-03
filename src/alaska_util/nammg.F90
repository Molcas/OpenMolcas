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

subroutine NAMmG( &
#                define _CALLING_
#                include "mem_interface.fh"
                )

use Definitions, only: iwp

implicit none
#define _USE_WP_
#include "mem_interface.fh"
integer(kind=iwp) :: iAng(4)
! Statement function
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

#include "macros.fh"
unused_var(lr)

iAng(1) = la
iAng(2) = lb
iAng(3) = 0
iAng(4) = 0
call MemRg1(iAng,nHer,Mem)
Mem = Mem+2+nElem(la)*nElem(lb)

return

end subroutine NAMmG
