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

use casvb_global, only: cvb, iorts, norb, nort, npr, nprorb, nvb, orbopt, orbs, sorbs, strucopt, sym, vbdet
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: dx(*)
integer(kind=iwp) :: ic
real(kind=wp), allocatable :: orbs1(:,:), cvb1(:)
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

if (orbopt) call touch_cvb('ORBS')
if (strucopt) call touch_cvb('CVB')
call make_cvb('WFN')
! TRY quantities only up2date if SVB/EVB ok:
if (up2date_cvb('SVBTRY')) call make_cvb('SVB')
if (up2date_cvb('EVBTRY')) call make_cvb('EVB')

ic = 1
call mma_allocate(orbs1,norb,norb,label='orbs1')
call mma_allocate(cvb1,nvb,label='cvb1')
call update2_cvb(orbs1,cvb1,orbs,cvb,sorbs,dx,ic,norb,nvb,nprorb,npr,orbopt,strucopt,sym,iorts,nort)
orbs(:,:) = orbs1(:,:)
cvb(1:nvb) = cvb1(:)
call str2vbc_cvb(cvb,vbdet)
call mma_deallocate(orbs1)
call mma_deallocate(cvb1)

return

end subroutine update_cvb
