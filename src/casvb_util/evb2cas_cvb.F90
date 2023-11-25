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

subroutine evb2cas_cvb(orbs,cvb,fx,ioptc,iter)

use casvb_global, only: civb1, civb2, civb3, civb4, civb5, norb, nvb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_IN_) :: orbs(norb,norb), cvb(nvb)
real(kind=wp), intent(out) :: fx
integer(kind=iwp), intent(out) :: ioptc, iter
integer(kind=iwp) :: idum
real(kind=wp) :: dx_amx, dxnrm
real(kind=wp), allocatable :: tmp(:)
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: tstfile_cvb ! ... Files/Hamiltonian available ...

if (tstfile_cvb(66000.2_wp)) then
  call mma_allocate(tmp,norb*norb+nvb,label='tmp')
  call rdlow_cvb(tmp,norb*norb+nvb,66000.2_wp,0)
  tmp(1:norb*norb) = pack(orbs(:,:),.true.)-tmp(1:norb*norb)
  tmp(norb*norb+1:) = cvb(:)-tmp(norb*norb+1:)
  dxnrm = dnrm2_(norb*norb+nvb,tmp,1)
  call findamx_cvb(tmp,norb*norb+nvb,dx_amx,idum)
  call mma_deallocate(tmp)
end if
call wrlow_cvb(orbs,norb*norb,66000.2_wp,0)
call wrlow_cvb(cvb,nvb,66000.2_wp,norb*norb)
call evb2cas2_cvb(orbs,cvb,ioptc,iter,fx,dxnrm,dx_amx,civb1,civb2,civb3,civb4,civb5)

return

end subroutine evb2cas_cvb
