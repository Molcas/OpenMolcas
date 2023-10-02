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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine update_cvb(dx)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: dx(*)
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: ic
real(kind=wp), allocatable :: orbs(:,:), cvb(:)
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

if (orbopt) call touch_cvb('ORBS')
if (strucopt) call touch_cvb('CVB')
call make_cvb('WFN')
! TRY quantities only up2date if SVB/EVB ok:
if (up2date_cvb('SVBTRY')) call make_cvb('SVB')
if (up2date_cvb('EVBTRY')) call make_cvb('EVB')

ic = 1
call mma_allocate(orbs,norb,norb,label='orbs')
call mma_allocate(cvb,nvb,label='nvb')
call update2_cvb(orbs,cvb,work(lv(1)),work(lv(2)),work(lw(2)),dx,ic,norb,nvb,nprorb,npr,orbopt,strucopt,sym,work(lp(6)), &
                 iwork(ls(11)),nort)
call fmove_cvb(orbs,work(lv(1)),norb*norb)
call fmove_cvb(cvb,work(lv(2)),nvb)
call str2vbc_cvb(work(lv(2)),work(lv(5)))
call mma_deallocate(orbs)
call mma_deallocate(cvb)

return

end subroutine update_cvb
