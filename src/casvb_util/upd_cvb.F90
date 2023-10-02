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

subroutine upd_cvb(dx,orbs,cvb)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: dx(*), orbs(norb,norb), cvb(nvb)
#include "WrkSpc.fh"
integer(kind=iwp) :: ic

if (orbopt) call touch_cvb('ORBSTRY')
if (strucopt) call touch_cvb('CVBTRY')
call make_cvb('WFNTRY')
ic = 2
call update2_cvb(orbs,cvb,work(lv(1)),work(lv(2)),work(lw(2)),dx,ic,norb,nvb,nprorb,npr,orbopt,strucopt,sym,work(lp(6)), &
                 iwork(ls(11)),nort)

return

end subroutine upd_cvb
