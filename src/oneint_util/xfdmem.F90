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

use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: iAngV(4), kab, lab, labcd, labMax, labMin, lc, lcd, lcdMax, lcdMin, ld, Mem1, Mem2, nFlop, nMem

lc = lr
ld = 0
nHer = (la+lb+lc+ld+2)/2
labMin = nTri3_Elem1(max(la,lb)-1)
labMax = nTri3_Elem1(la+lb)-1
lcdMin = nTri3_Elem1(lr-1)
lcdMax = nTri3_Elem1(lr)-1
lab = (labMax-labMin+1)
kab = nTri_Elem1(la)*nTri_Elem1(lb)
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
