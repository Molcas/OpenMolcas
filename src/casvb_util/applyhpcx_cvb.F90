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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: civec(*), c_daxpy
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: ibase, icivec, isyml, isymmx, nci
real(kind=wp) :: cnrm
real(kind=wp), allocatable :: cim(:), cim2(:)
real(kind=wp), parameter :: thr2 = 1.0e-20_wp
integer(kind=iwp), external :: mstackr_cvb
real(kind=wp), external :: ddot_

icivec = nint(civec(1))
n_applyh = n_applyh+1
call setcnt2_cvb(icivec,0)
if (iform_ci(icivec) /= 0) then
  write(u6,*) ' Unsupported format in APPLYH :',iform_ci(icivec)
  call abend_cvb()
end if

! (NIRREP may be altered in loop)
isymmx = nirrep
do isyml=1,isymmx
  nci = ncivb(isyml)
  ibase = mstackr_cvb(0)
  call mma_allocate(cim,nci,label='cim')
  cim(:) = Zero
  ibasemx = max(ibasemx,mstackr_cvb(0))
  call vb2mol_cvb(work(iaddr_ci(icivec)),cim,isyml)

  ! If only one irrep present keep down memory requirements:
  if ((isymmx > 1) .and. (nci /= ndet)) then
    call mma_allocate(cim2,nci,label='cim2')
    cim2(:) = Zero
    ibasemx = max(ibasemx,mstackr_cvb(0))
    cnrm = ddot_(nci,cim,1,cim,1)
    ! If anything there, apply Hamiltonian to vector of this symmetry:
    if (cnrm > thr2) call sigmadet_cvb(cim,cim2,isyml,nci)
    if (c_daxpy /= Zero) call daxpy_(nci,c_daxpy,cim,1,cim2,1)
    call mol2vb_cvb(work(iaddr_ci(icivec)),cim2,isyml)
    call mma_deallocate(cim2)
  else
    call fzero(work(iaddr_ci(icivec)),nci)
    cnrm = ddot_(nci,cim,1,cim,1)
    ! If anything there, apply Hamiltonian to vector of this symmetry:
    if (cnrm > thr2) then
      call fzero(work(iaddr_ci(icivec)),nci)
      call sigmadet_cvb(cim,work(iaddr_ci(icivec)),isyml,nci)
    end if
    if (c_daxpy /= Zero) call daxpy_(nci,c_daxpy,cim,1,work(iaddr_ci(icivec)),1)
    call fmove_cvb(work(iaddr_ci(icivec)),cim,nci)
    call mol2vb_cvb(work(iaddr_ci(icivec)),cim,isyml)
  end if
  call mma_deallocate(cim)
  call mfreer_cvb(ibase)
end do

return

end subroutine applyhpcx_cvb
