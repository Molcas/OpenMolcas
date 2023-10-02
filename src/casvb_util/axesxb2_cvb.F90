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

subroutine axesxb2_cvb(asonc,ddres2upd,vec,resthr_inp,ioptc,iter,fx_exp,c,axc,sxc,share,res,rhs,ap,rhsp,solp,solp_res)

use casvb_global, only: corenrg, ifollow, ipdd, isaddledd, maxd, mxit, nfrdim, nortiterdd, nparm, nvguess, nvrestart, orththrdd, &
                        resthrdd
use casvb_interfaces, only: ddasonc_sub, ddres_sub, ddres2upd_sub, ddrestart_sub, ddsol_sub
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
procedure(ddasonc_sub) :: asonc
procedure(ddres2upd_sub) :: ddres2upd
real(kind=wp) :: vec(nparm), resthr_inp, fx_exp, c(nparm,maxd), axc(nparm,maxd), sxc(nparm,maxd), res(nparm), rhs(*), &
                 ap(maxd,maxd), rhsp(maxd), solp(maxd), solp_res(maxd)
integer(kind=iwp) :: ioptc, iter
logical(kind=iwp) :: share
real(kind=wp) :: resthr_use
procedure(ddres_sub) :: axesxbres_cvb
procedure(ddrestart_sub) :: ddrestart_cvb
procedure(ddsol_sub) :: axexbsol_cvb

! If RESTHR_INP unset use default:
if (resthr_inp /= Zero) then
  resthr_use = resthr_inp
else
  resthr_use = resthrdd
end if
!------------------------------------------
! ASonC             ASonC10_CVB
! AxExbSOL_CVB      DDSOL10_CVB
! AxESxbRES_CVB     DDRES10_CVB
! DDRES2UPD         DDRES2UPD10_CVB
! DDRESTART_CVB     DDRESTART_CVB
!------------------------------------------
call dirdiag_cvb(asonc,axexbsol_cvb,axesxbres_cvb,ddres2upd,ddrestart_cvb,c,axc,sxc,share,vec,res,rhs,ap,rhsp,solp,solp_res, &
                 .false.,.true.,.true.,maxd,nparm,nfrdim,nvguess,nvrestart,isaddledd,ifollow,mxit,resthr_use,orththrdd,nortiterdd, &
                 corenrg,ioptc,iter,fx_exp,ipdd)

return

end subroutine axesxb2_cvb
