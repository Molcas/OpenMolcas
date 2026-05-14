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

subroutine applyhpcx_cvb(civec,c_daxpy)
! Exact copy if applyh except for c_daxpy in arg list.

use casvb_global, only: icnt_ci, iform_ci, n_applyh, ncivb, ndet, nirrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: civec(0:ndet)
real(kind=wp), intent(in) :: c_daxpy
integer(kind=iwp) :: icivec, isyml, isymmx, nci
real(kind=wp) :: cnrm
real(kind=wp), allocatable :: cim(:), cim2(:)
real(kind=wp), parameter :: thr2 = 1.0e-20_wp
real(kind=wp), external :: ddot_

n_applyh = n_applyh+1
icivec = nint(civec(0))
icnt_ci(icivec) = 0
if (iform_ci(icivec) /= 0) then
  write(u6,*) ' Unsupported format in APPLYH :',iform_ci(icivec)
  call abend_cvb()
end if

! (NIRREP may be altered in loop)
isymmx = nirrep
do isyml=1,isymmx
  nci = ncivb(isyml)
  call mma_allocate(cim,nci,label='cim')
  cim(:) = Zero
  call vb2mol_cvb(civec(1:),cim,isyml)

  ! If only one irrep present keep down memory requirements:
  if ((isymmx > 1) .and. (nci /= ndet)) then
    call mma_allocate(cim2,nci,label='cim2')
    cim2(:) = Zero
    cnrm = ddot_(nci,cim,1,cim,1)
    ! If anything there, apply Hamiltonian to vector of this symmetry:
    if (cnrm > thr2) call sigmadet_cvb(cim,cim2,isyml,nci)
    if (c_daxpy /= Zero) cim2(:) = cim2(:)+c_daxpy*cim(:)
    call mol2vb_cvb(civec(1:),cim2,isyml)
    call mma_deallocate(cim2)
  else
    civec(1:nci) = Zero
    cnrm = ddot_(nci,cim,1,cim,1)
    ! If anything there, apply Hamiltonian to vector of this symmetry:
    if (cnrm > thr2) call sigmadet_cvb(cim,civec(1:),isyml,nci)
    if (c_daxpy /= Zero) civec(1:nci) = civec(1:nci)+c_daxpy*cim(:)
    cim(:) = civec(1:nci)
    call mol2vb_cvb(civec(1:),cim,isyml)
  end if
  call mma_deallocate(cim)
end do

return

end subroutine applyhpcx_cvb
