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

subroutine XFdMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

#include "mem_interface.fh"
integer iAngV(4)
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

lc = lr
ld = 0
nHer = (la+lb+lc+ld+2)/2
labMin = nabSz(max(la,lb)-1)+1
labMax = nabSz(la+lb)
lcdMin = nabSz(lr-1)+1
lcdMax = nabSz(lr)
lab = (labMax-labMin+1)
kab = nElem(la)*nElem(lb)
lcd = (lcdMax-lcdMin+1)
labcd = lab*lcd

call mHRR(la,lb,nFlop,nMem)
Mem1 = max(lcd*nMem,labcd)

iAngV(1) = la
iAngV(2) = lb
iAngV(3) = lc
iAngV(4) = ld
call MemRys(iAngV,Mem2)
Mem2 = max(Mem2,kab*lcd)

Mem = Mem1+Mem2

return

end subroutine XFdMem
