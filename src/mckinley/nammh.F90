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

subroutine NAMmh(nRys,MmNAG,la,lb,lr)

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nRys, MmNAG, la, lb, lr
integer iAng(4)

iAng(1) = la
iAng(2) = lb
iAng(3) = 0
iAng(4) = 0
call MemRg2(iAng,nRys,MmNAG,2)
MmNAG = MmNAG+2+nTri_Elem1(la)*nTri_Elem1(lb) ! Alpha beta & DAO

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(lr)

end subroutine NAMmh
