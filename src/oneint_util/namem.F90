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

subroutine NAMem( &
#                define _CALLING_
#                include "mem_interface.fh"
                )

use Basis_Info, only: mGaussian_Type, Nuclear_Model
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: iAngV(4), labcd, MemNA2, nFlop, nMem

!                                                                      *
!***********************************************************************
!                                                                      *
call mHrr(la,lb,nFlop,nMem)

iAngV(1) = la
iAngV(2) = lb
iAngV(3) = lr
iAngV(4) = 0
call MemRys(iAngV,Mem)
nHer = (la+lb+lr+2)/2
if (Nuclear_Model == mGaussian_Type) then

  labcd = (la+1)*(la+2)/2*(lb+1)*(lb+2)/2

  iAngV(3) = lr+2
  call MemRys(iAngV,MemNA2)
  Mem = max(Mem,MemNA2)
  nHer = (la+lb+lr+4)/2
  Mem = Mem+labcd
end if

Mem = max(nMem,Mem)

return

end subroutine NAMem
