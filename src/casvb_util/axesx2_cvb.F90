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

subroutine axesx2_cvb(asonc,ddres2upd,vec,resthr_inp,ioptc,iter,fx_exp,c,axc,sxc,share,res,ap,solp,solp_res)

use casvb_global, only: corenrg, ifollow, ipdd, isaddledd, maxd, mxit, nfrdim, nortiterdd, nparm, nvguess, nvrestart, orththrdd, &
                        resthrdd
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
external :: asonc, ddres2upd
real(kind=wp) :: vec(nparm), resthr_inp, fx_exp, c(nparm,maxd), axc(nparm,maxd), sxc(nparm,maxd), res(nparm), ap(maxd,maxd), &
                 solp(maxd), solp_res(maxd)
integer(kind=iwp) :: ioptc, iter
logical(kind=iwp) :: share
real(kind=wp) :: dum(1), resthr_use
external :: axesxres_cvb, axexsol_cvb, ddrestart_cvb

! If RESTHR_INP unset use default:
if (resthr_inp /= Zero) then
  resthr_use = resthr_inp
else
  resthr_use = resthrdd
end if

call dirdiag_cvb(asonc,axexsol_cvb,axesxres_cvb,ddres2upd,ddrestart_cvb,c,axc,sxc,share,vec,res,dum,ap,dum,solp,solp_res,.false., &
                 .true.,.false.,maxd,nparm,nfrdim,nvguess,nvrestart,isaddledd,ifollow,mxit,resthr_use,orththrdd,nortiterdd, &
                 corenrg,ioptc,iter,fx_exp,ipdd)

return

end subroutine axesx2_cvb
