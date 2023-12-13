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

subroutine upd_cvb(dx,orbs1,cvb1)

use casvb_global, only: cvb, iorts, norb, nort, npr, nprorb, nvb, orbopt, orbs, sorbs, strucopt, sym
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: dx(*)
real(kind=wp), intent(out) :: orbs1(norb,norb), cvb1(nvb)
integer(kind=iwp) :: ic

if (orbopt) call touch_cvb('ORBSTRY')
if (strucopt) call touch_cvb('CVBTRY')
call make_cvb('WFNTRY')
ic = 2
call update2_cvb(orbs1,cvb1,orbs,cvb,sorbs,dx,ic,norb,nvb,nprorb,npr,orbopt,strucopt,sym,iorts,nort)

return

end subroutine upd_cvb
