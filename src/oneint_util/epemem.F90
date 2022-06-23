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

subroutine EPEMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: iAngV(4), MemHrr, nFlop

#include "macros.fh"
unused_var(lr)

call mHrr(la,lb,nFlop,MemHrr)

nHer = (la+lb+2)/2
iAngV(1) = la
iAngV(2) = lb
iAngV(3) = 0
iAngV(4) = 0
call MemRys(iAngV,Mem)

Mem = max(Mem,MemHrr)

return

end subroutine EPEMem
