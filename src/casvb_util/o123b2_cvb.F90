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

use casvb_global, only: cnrm, isaddle, maxize, safety, signtol
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nparm
real(kind=wp) :: dx(nparm), eigvec(nparm,nparm), eigval(nparm), dxp(nparm), gradp(nparm), wrk(nparm), dxnrm
integer(kind=iwp) :: nnegeig, nposeig
real(kind=wp) :: alfastart, eig, eigmn, eigmx, safety_use
real(kind=wp), external :: dnrm2_

if (maxize) then
  nposeig = min(isaddle,nparm)
else
  nposeig = max(nparm-isaddle,nparm)
end if
nnegeig = nparm-nposeig

eigmx = -One
eigmn = One
if (nnegeig > 0) eigmx = eigval(nnegeig)
if (nposeig > 0) eigmn = eigval(nnegeig+1)
safety_use = safety
!do
if ((eigmx < -signtol) .and. (eigmn > signtol)) then
  alfastart = Zero
else
  alfastart = max(eigmx,-eigmn,Zero)+safety_use
end if
call getdxp_cvb(dxp,gradp,eigval,nnegeig,nparm,alfastart)
cnrm = dnrm2_(nparm,dxp,1)
!  Increased level shift not necessary when full Hessian is calculated:
!  if (alfastart == Zero) exit
!  gnrm = dnrm2_(nparm,gradp,1)
!  if ((cnrm <= 1.0e-15_wp) .or. (gnrm <= 1.0e-15_wp) .or. (safety_use == 1.0e-4_wp)) exit
!  ovr_dx_grad = ddot_(nparm,dxp,1,gradp,1)/(cnrm*gnrm)
!  if (ovr_dx_grad >= 0.3_wp) exit
!  safety_use = 1.0e-4_wp
!end do

call makedx_cvb(dx,nparm,0,eigvec,eigval,dxp,gradp,wrk,.false.,.false.,nposeig,.false.,.false.,nnegeig,.false.,alfastart,eig)
dxnrm = dnrm2_(nparm,dx,1)

return

end subroutine o123b2_cvb
