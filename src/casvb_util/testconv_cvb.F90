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

subroutine testconv_cvb(fx,npr,dx,w2,exp_tc,close2conv,converged,wrongstat)

use casvb_global, only: fxbest, ip, isaddleo, maxize
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: npr
real(kind=wp), intent(in) :: fx, dx(npr), w2(npr), exp_tc
logical(kind=iwp), intent(inout) :: close2conv
logical(kind=iwp), intent(out) :: converged, wrongstat
integer(kind=iwp) :: nnegeig, nposeig
real(kind=wp) :: act, eigmn, eigmna, eigmx, zz

if (maxize) then
  nposeig = min(isaddleo,npr)
else
  nposeig = max(npr-isaddleo,npr)
end if
nnegeig = npr-nposeig

! EIGMX is maximum of NNEGEIG first Hessian eigenvalues (which should
! all be negative) EIGMN the minimum of NPOSEIG last eigenvalues:
eigmx = -One
eigmn = One
eigmna = One

call zz_cvb(act,zz,fx,fxbest,exp_tc,ip)
fxbest = fx

! << Test for convergence or near-convergence >>
call testconv2_cvb(close2conv,converged,wrongstat,act,zz,dx,w2,npr,eigmn,eigmx,eigmna,nposeig,nnegeig)

return

end subroutine testconv_cvb
