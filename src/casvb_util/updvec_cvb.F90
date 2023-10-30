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

subroutine updvec_cvb(upd,iorb,jorb,niprev,iprev,orbs,north,corth)
! Find update for IORB as projection of JORB on allowed space

use casvb_global, only: niorth, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: upd(norb)
integer(kind=iwp), intent(in) :: iorb, jorb, niprev, iprev(niprev), north(norb)
real(kind=wp), intent(in) :: orbs(norb,norb), corth(norb,niorth)
integer(kind=iwp) :: i, ncon, noffort
real(kind=wp) :: dum(1)
real(kind=wp), allocatable :: tmp(:,:)

call mma_allocate(tmp,norb,norb,label='tmp')
noffort = sum(north(1:iorb-1))
! Collect all constraints and find span:
call span0_cvb(norb,norb)
if (north(iorb) > 0) call span1_cvb(corth(:,1+noffort),north(iorb),dum,norb,0)
do i=1,niprev
  call span1_cvb(orbs(:,iprev(i)),1,dum,norb,0)
end do
call span1_cvb(orbs(:,iorb),1,dum,norb,0)
call span2_cvb(tmp,ncon,dum,norb,0)

! Orthogonalise update to all remaining constraints
upd(:) = orbs(:,jorb)
call schmidtd_cvb(tmp,ncon,upd,1,dum,norb,0)
call mma_deallocate(tmp)

return

end subroutine updvec_cvb
