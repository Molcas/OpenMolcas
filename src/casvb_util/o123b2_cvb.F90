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

subroutine o123b2_cvb(nparm,dx,eigvec,eigval,dxp,gradp,wrk,dxnrm)

implicit real*8(a-h,o-z)
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"
dimension dx(nparm)
dimension eigvec(nparm,nparm), eigval(nparm)
dimension dxp(nparm), gradp(nparm), wrk(nparm)
save zero, one
data zero/0.d0/,one/1d0/

if (maxize) then
  nposeig = min(isaddle,nparm)
else
  nposeig = max(nparm-isaddle,nparm)
end if
nnegeig = nparm-nposeig

eigmx = -one
eigmn = one
if (nnegeig > 0) eigmx = eigval(nnegeig)
if (nposeig > 0) eigmn = eigval(nnegeig+1)
safety_use = safety
!do
if ((eigmx < -signtol) .and. (eigmn > signtol)) then
  alfastart = zero
else
  alfastart = max(eigmx,-eigmn,zero)+safety_use
end if
call getdxp_cvb(dxp,gradp,eigval,nnegeig,nparm,alfastart)
cnrm = dnrm2_(nparm,dxp,1)
!  Increased level shift not necessary when full Hessian is calculated:
!  if (alfastart == zero) exit
!  gnrm = dnrm2_(nparm,gradp,1)
!  if ((cnrm <= 1d-15) .or. (gnrm <= 1d-15) .or. (safety_use == 1d-4)) exit
!  ovr_dx_grad = ddot_(nparm,dxp,1,gradp,1)/(cnrm*gnrm)
!  if (ovr_dx_grad >= .3d0) exit
!  safety_use = 1d-4
!end do

call makedx_cvb(dx,nparm,0,eigvec,eigval,dxp,gradp,wrk,.false.,.false.,nposeig,.false.,.false.,nnegeig,.false.,alfastart,eig)
dxnrm = dnrm2_(nparm,dx,1)

return

end subroutine o123b2_cvb
