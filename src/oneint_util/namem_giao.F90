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

subroutine NAMem_GIAO( &
#                     define _CALLING_
#                     include "mem_interface.fh"
                     )

use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: iAngV(4), kab, lab, labcd_EF, labMax, labMin, lc, lcd_EF, lcd_NA, lcdMax_EF, lcdMax_NA, lcdMin_EF, lcdMin_NA, &
                     ld, Mem0, Mem1, Mem2, Mem2_EF, Mem2_NA, nFlop, nMem

lc = lr
ld = 0
nHer = (la+lb+lc+ld+2)/2

labMin = nTri3_Elem1(max(la,lb)-1)
labMax = nTri3_Elem1(la+lb)-1
lab = (labMax-labMin+1)
kab = nTri_Elem1(la)*nTri_Elem1(lb)

lcdMin_EF = nTri3_Elem1(lr-1)
lcdMax_EF = nTri3_Elem1(lr)-1
lcd_EF = (lcdMax_EF-lcdMin_EF+1)
labcd_EF = lab*lcd_EF
Mem0 = labcd_EF

lcdMin_NA = nTri3_Elem1(lr-2)
lcdMax_NA = nTri3_Elem1(lr-1)-1
lcd_NA = (lcdMax_NA-lcdMin_NA+1)

call mHRR(la,lb,nFlop,nMem)
Mem1 = max(lcd_EF,lcd_NA)*nMem

iAngV(1) = la
iAngV(2) = lb
iAngV(3) = lc
iAngV(4) = ld
call MemRys(iAngV,Mem2_EF)
iAngV(3) = 0
call MemRys(iAngV,Mem2_NA)
Mem2 = max(Mem2_EF,Mem2_NA,kab*max(lcd_EF,lcd_NA))

Mem = Mem0+Mem1+Mem2

return

end subroutine NAMem_GIAO
