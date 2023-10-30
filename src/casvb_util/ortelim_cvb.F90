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

subroutine ortelim_cvb(trprm,iorts,irots,sorbs,nc,npr1,norbprm,nrem)

use casvb_global, only: ndrot, norb, nort
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iorts(2,nort), irots(2,ndrot), nc, npr1, norbprm
real(kind=wp), intent(inout) :: trprm(npr1,npr1)
real(kind=wp), intent(in) :: sorbs(norb,norb)
integer(kind=iwp), intent(out) :: nrem
integer(kind=iwp) :: iorb, iort, irot, jorb, ki, kj, korb
real(kind=wp) :: dum(1)
real(kind=wp), allocatable :: tmp(:,:)

call mma_allocate(tmp,norbprm,max(nc+nort+ndrot,norbprm),label='tmp')
tmp(:,:) = Zero
tmp(:,1:nc) = trprm(1:norbprm,1:nc)
do iort=1,nort
  iorb = iorts(1,iort)
  jorb = iorts(2,iort)
  do korb=1,norb
    ki = korb+(iorb-1)*(norb-1)
    if (korb > iorb) ki = ki-1
    kj = korb+(jorb-1)*(norb-1)
    if (korb > jorb) kj = kj-1
    if (korb /= iorb) tmp(ki,iort+nc) = sorbs(korb,jorb)
    if (korb /= jorb) tmp(kj,iort+nc) = sorbs(korb,iorb)
  end do
end do
do irot=1,ndrot
  iorb = irots(1,irot)
  jorb = irots(2,irot)
  do korb=1,norb
    ki = korb+(iorb-1)*(norb-1)
    if (korb > iorb) ki = ki-1
    kj = korb+(jorb-1)*(norb-1)
    if (korb > jorb) kj = kj-1
    if (korb /= iorb) tmp(ki,irot+nc+nort) = sorbs(korb,jorb)
    if (korb /= jorb) tmp(kj,irot+nc+nort) = -sorbs(korb,iorb)
  end do
end do
call span_cvb(tmp,nc+nort+ndrot,nrem,dum,norbprm,0)
call compl_cvb(tmp,nrem,norbprm)

trprm(:,:) = Zero
trprm(1:norbprm,1:norbprm) = tmp(1:norbprm,1:norbprm)

call mma_deallocate(tmp)

return

end subroutine ortelim_cvb
