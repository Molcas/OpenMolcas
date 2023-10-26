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

subroutine o12sb_cvb( &
#                    define _CALLING_
#                    include "optb_interface.fh"
                    )

use casvb_global, only: cvb, expct, fxbest, have_solved_it, hh, ip, nfrorb, nvb, odx, orbs, scalesmall, strucopt
use casvb_interfaces, only: ddasonc_sub, ddres2upd_sub
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
#include "optb_interface.fh"
integer(kind=iwp) :: ioptc, ipu, iter
real(kind=wp) :: cnrm2, fac, fx_exp, resthr_old = -One, resthr_use
logical(kind=iwp) :: skip
real(kind=wp), external :: ddot_, dnrm2_
procedure(ddasonc_sub) :: asonc12s_cvb
procedure(ddres2upd_sub) :: ddres2upd10_cvb

if (.not. close2conv) then
  resthr_use = 1.0e-5_wp
else
  resthr_use = 0.05_wp*grdnrm
  resthr_use = min(1.0e-5_wp,resthr_use)
  resthr_use = max(1.0e-9_wp,resthr_use)
end if
skip = ((resthr_use == resthr_old) .and. have_solved_it)
resthr_old = resthr_use

if (.not. skip) then
  call makegjorbs_cvb(orbs)

  call axb_cvb(asonc12s_cvb,ddres2upd10_cvb,odx,resthr_use,ioptc,iter,fx_exp)
  expct = fx_exp-fxbest
  have_solved_it = .true.

  if (ip >= 2) write(u6,'(a,i4)') ' Number of iterations for direct diagonalization :',iter

  if (strucopt) then
    cnrm2 = ddot_(nvb,cvb,1,odx(nfrorb+1:),1)
    ! "Orthogonalize" on CVB to get smallest possible update norm:
    odx(nfrorb+1:nfrorb+nvb) = odx(nfrorb+1:nfrorb+nvb)-cnrm2*cvb(1:nvb)
    ! Scale variables according to overlap with CVB:
    odx(1:nparm) = odx(1:nparm)/cnrm2
  else
    ! We are doing "Augmented" calc:
    fac = One/odx(1)
    ! Scale variables according to overlap with CVB:
    odx(1:nparm-1) = fac*odx(2:nparm)
  end if
end if

dxnrm = dnrm2_(nparm,odx,1)
if (.not. close2conv) then
  ipu = 1
else
  ipu = 2
end if
if ((dxnrm > hh) .or. scalesmall(ipu)) then
  odx(1:nparm) = hh/dxnrm*odx(1:nparm)
  dxnrm = hh
end if

return

end subroutine o12sb_cvb
