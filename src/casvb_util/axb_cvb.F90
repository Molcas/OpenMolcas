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

subroutine axb_cvb(asonc,ddres2upd,vec,resthr_inp,ioptc,iter,fx_exp)

use casvb_global, only: c, corenrg, ifollow, ipdd, isaddledd, maxd, mxit, nfrdim, nortiterdd, nparm, nvguess, nvrestart, &
                        orththrdd, res, resthrdd, rhs, rhsp, solp, solp_res, sxc
use casvb_interfaces, only: ddasonc_sub, ddres_sub, ddres2upd_sub, ddsol_sub
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
procedure(ddasonc_sub) :: asonc
procedure(ddres2upd_sub) :: ddres2upd
real(kind=wp), intent(out) :: vec(nparm), fx_exp
real(kind=wp), intent(in) :: resthr_inp
integer(kind=iwp), intent(out) :: ioptc, iter
real(kind=wp) :: dum(1), resthr_use
procedure(ddres_sub) :: axbres_cvb
procedure(ddsol_sub) :: axbsol_cvb

! If RESTHR_INP unset use default:
if (resthr_inp /= Zero) then
  resthr_use = resthr_inp
else
  resthr_use = resthrdd
end if
call dirdiag_cvb(asonc,axbsol_cvb,axbres_cvb,ddres2upd,c,dum,sxc,.false.,vec,res,rhs,dum,rhsp,solp,solp_res,.false.,.false., &
                 .true.,maxd,nparm,nfrdim,nvguess,nvrestart,isaddledd,ifollow,mxit,resthr_use,orththrdd,nortiterdd,corenrg,ioptc, &
                 iter,fx_exp,ipdd)

return

end subroutine axb_cvb
