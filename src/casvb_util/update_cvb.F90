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

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: dx(*)
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, i2, i3, ic
integer(kind=iwp), external :: mstackr_cvb
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

if (orbopt) call touch_cvb('ORBS')
if (strucopt) call touch_cvb('CVB')
call make_cvb('WFN')
! TRY quantities only up2date if SVB/EVB ok:
if (up2date_cvb('SVBTRY')) call make_cvb('SVB')
if (up2date_cvb('EVBTRY')) call make_cvb('EVB')

ic = 1
i1 = mstackr_cvb(norb*norb)
i2 = mstackr_cvb(nvb)
i3 = mstackr_cvb(norb*norb)
call update2_cvb(work(i1),work(i2),work(lv(1)),work(lv(2)),work(lw(2)),dx,ic,norb,nvb,nprorb,npr,orbopt,strucopt,sym,work(lp(6)), &
                 iwork(ls(11)),nort,work(i3))
call fmove_cvb(work(i1),work(lv(1)),norb*norb)
call fmove_cvb(work(i2),work(lv(2)),nvb)
call str2vbc_cvb(work(lv(2)),work(lv(5)))
call mfreer_cvb(i1)

return

end subroutine update_cvb
