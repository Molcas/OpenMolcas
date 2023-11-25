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

subroutine axexbsol_cvb( &
#                       define _CALLING_
#                       include "ddsol_interface.fh"
                       )
! Solve linear equation in Davidson subspace.

use casvb_global, only: cnrm, corenrg, ifollow, ipdd, nroot, safety, signtol
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "ddsol_interface.fh"
integer(kind=iwp) :: i, ioptc, nnegeig, nposeig
real(kind=wp) :: alfastart, eigmn, eigmx, gnrm, ovr_dx_grad, safety_use
real(kind=wp), allocatable :: dxp(:), eigval(:), eigvec(:,:), gradp(:), w2(:)
real(kind=wp), external :: ddot_, dnrm2_

#include "macros.fh"
unused_var(nfrdim)

call mma_allocate(eigval,itdav,label='eigval')
call mma_allocate(eigvec,itdav,itdav,label='eigvec')
eigvec(:,:) = ap(1:itdav,1:itdav)

if (ipdd >= 3) then
  write(u6,*) ' AP matrix :'
  do i=1,itdav
    eigval(i) = eigvec(i,i)
    eigvec(i,i) = eigvec(i,i)+corenrg
  end do
  call mxprintd_cvb(eigvec,itdav,itdav,0)
  do i=1,itdav
    eigvec(i,i) = eigval(i)
  end do
  write(u6,*) ' RHSP vector :'
  call mxprint_cvb(rhsp,1,itdav,0)
end if

call mxdiag_cvb(eigvec,eigval,itdav)

if (ipdd >= 2) then
  write(u6,'(a)') ' Eigenvalues :'
  eigval(1:itdav) = eigval(1:itdav)+corenrg
  call vecprint_cvb(eigval,itdav)
  eigval(1:itdav) = eigval(1:itdav)-corenrg
end if

call mma_allocate(gradp,itdav,label='gradp')
call mxatb_cvb(rhsp,eigvec,1,itdav,itdav,gradp)
if (ifollow == 1) then
  nposeig = nroot-1
  nnegeig = itdav-nposeig
else if (ifollow == 2) then
  nnegeig = nroot-1
  nposeig = itdav-nnegeig
else
  write(u6,*) ' Error in IFOLLOW with direct Fletcher!',ifollow
  call abend_cvb()
  nnegeig = 0
  nposeig = 0
end if

call mma_allocate(dxp,itdav,label='dxp')
eigmx = -One
eigmn = One
if (nnegeig > 0) eigmx = eigval(nnegeig)
if (nposeig > 0) eigmn = eigval(nnegeig+1)
safety_use = safety
do
  if ((eigmx < -signtol) .and. (eigmn > signtol)) then
    alfastart = Zero
  else
    alfastart = max(eigmx,-eigmn,Zero)+safety_use
  end if
  call getdxp_cvb(dxp,gradp,eigval,nnegeig,itdav,alfastart)
  cnrm = dnrm2_(itdav,dxp,1)
  if (alfastart == Zero) exit
  gnrm = dnrm2_(itdav,gradp,1)
  if ((cnrm <= 1.0e-15_wp) .or. (gnrm <= 1.0e-15_wp) .or. (safety_use == 1.0e-4_wp)) exit
  ovr_dx_grad = ddot_(itdav,dxp,1,gradp,1)/(cnrm*gnrm)
  if (ovr_dx_grad >= 0.3_wp) exit
  safety_use = 1.0e-4_wp
end do

call mma_allocate(w2,itdav,label='w2')
ioptc = 0
call makedx_cvb(solp,itdav,ioptc,eigvec,eigval,dxp,gradp,w2,.false.,nposeig,.false.,.false.,nnegeig,.false.,alfastart,eig)
call mma_deallocate(eigval)
call mma_deallocate(eigvec)
call mma_deallocate(gradp)
call mma_deallocate(dxp)
call mma_deallocate(w2)

eig_res = eig
solp_res(:) = solp(:)
if (ipdd >= 2) then
  write(u6,'(a,f15.8)') ' Eigenvalue :',eig
  write(u6,'(a)') ' Solution vector :'
  call vecprint_cvb(solp,itdav)
end if

return

end subroutine axexbsol_cvb
