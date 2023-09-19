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

subroutine axexbsol2_cvb(ap,rhsp,itdav,maxdav,nfrdim1,solp,solp_res,eig,eig_res,eigval,eigvec,dxp,gradp,w2)
! Solve linear equation in Davidson subspace.

implicit real*8(a-h,o-z)
#include "direct_cvb.fh"
#include "locopt1_cvb.fh"
#include "tune_cvb.fh"
dimension ap(maxdav,maxdav), rhsp(maxdav)
dimension solp(maxdav), solp_res(maxdav)
dimension eigval(itdav), eigvec(itdav,itdav)
dimension dxp(itdav), gradp(itdav), w2(itdav)
save zero, one
data zero/0d0/,one/1d0/

do it=1,itdav
  call fmove_cvb(ap(1,it),eigvec(1,it),itdav)
end do

if (ip >= 3) then
  write(6,*) ' AP matrix :'
  do i=1,itdav
    eigval(i) = eigvec(i,i)
    eigvec(i,i) = eigvec(i,i)+corenrg
  end do
  call mxprintd_cvb(eigvec,itdav,itdav,0)
  do i=1,itdav
    eigvec(i,i) = eigval(i)
  end do
  write(6,*) ' RHSP vector :'
  call mxprint_cvb(rhsp,1,itdav,0)
end if

call mxdiag_cvb(eigvec,eigval,itdav)

if (ip >= 2) then
  write(6,'(a)') ' Eigenvalues :'
  do i=1,itdav
    eigval(i) = eigval(i)+corenrg
  end do
  call vecprint_cvb(eigval,itdav)
  do i=1,itdav
    eigval(i) = eigval(i)-corenrg
  end do
end if

call mxatb_cvb(rhsp,eigvec,1,itdav,itdav,gradp)
if (ifollow == 1) then
  nposeig = nroot-1
  nnegeig = itdav-nposeig
else if (ifollow == 2) then
  nnegeig = nroot-1
  nposeig = itdav-nnegeig
else
  write(6,*) ' Error in IFOLLOW with direct Fletcher!',ifollow
  call abend_cvb()
  nnegeig = 0
  nposeig = 0
end if

eigmx = -one
eigmn = one
if (nnegeig > 0) eigmx = eigval(nnegeig)
if (nposeig > 0) eigmn = eigval(nnegeig+1)
safety_use = safety
do
  if ((eigmx < -signtol) .and. (eigmn > signtol)) then
    alfastart = zero
  else
    alfastart = max(eigmx,-eigmn,zero)+safety_use
  end if
  call getdxp_cvb(dxp,gradp,eigval,nnegeig,itdav,alfastart)
  cnrm = dnrm2_(itdav,dxp,1)
  if (alfastart == zero) exit
  gnrm = dnrm2_(itdav,gradp,1)
  if ((cnrm <= 1d-15) .or. (gnrm <= 1d-15) .or. (safety_use == 1d-4)) exit
  ovr_dx_grad = ddot_(itdav,dxp,1,gradp,1)/(cnrm*gnrm)
  if (ovr_dx_grad >= .3d0) exit
  safety_use = 1d-4
end do

call makedx_cvb(solp,itdav,0,eigvec,eigval,dxp,gradp,w2,.false.,.false.,nposeig,.false.,.false.,nnegeig,.false.,alfastart,eig)

eig_res = eig
call fmove_cvb(solp,solp_res,itdav)
if (ip >= 2) then
  write(6,'(a,f15.8)') ' Eigenvalue :',eig
  write(6,'(a)') ' Solution vector :'
  call vecprint_cvb(solp,itdav)
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nfrdim1)

end subroutine axexbsol2_cvb
