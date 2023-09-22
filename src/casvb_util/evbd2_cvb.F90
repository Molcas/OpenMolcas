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

subroutine evbd2_cvb(orbs,cvb,fx,ioptc,iter,gjorb,gjorb2,gjorb3,c,axc,sxc,res,hp,solp,solp_res)

use casvb_global, only: follow, have_solved_it, nortiter, orththr, resthr
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
#include "optze_cvb.fh"
real(kind=wp) :: orbs(norb,norb), cvb(nvb), fx, gjorb(*), gjorb2(*), gjorb3(*), c(nvb,maxdav), axc(nvb,maxdav), sxc(nvb,maxdav), &
                 res(nvb), hp(maxdav,maxdav), solp(maxdav), solp_res(maxdav)
integer(kind=iwp) :: ioptc, iter
#include "print_cvb.fh"
integer(kind=iwp) :: ifollow, nvguess, nvrestart
real(kind=wp) :: dum(max(nvb,maxdav))
external :: asonc_cvb, ddres2upd10_cvb, ddres7_cvb, ddrestart_cvb, ddsol7_cvb

call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)

ioptc = 1
nvguess = 1
nvrestart = 0
call fmove_cvb(cvb,c,nvb)
if (.not. follow) then
  ifollow = 2
else
  ifollow = 4
end if

call ddinit7_cvb(ifollow,isaddle,ip(3))
call ddres2updinit_cvb(0)
call dirdiag_cvb(asonc_cvb,ddsol7_cvb,ddres7_cvb,ddres2upd10_cvb,ddrestart_cvb,c,axc,sxc,.false.,cvb,res,dum,hp,dum,solp,solp_res, &
                 .false.,.true.,.false.,maxdav,nvb,nvb,nvguess,nvrestart,isaddle,ifollow,mxiter,resthr,orththr,nortiter,corenrg, &
                 ioptc,iter,fx,ip(3))
have_solved_it = .true.

ovraa = one
evb = fx

return

end subroutine evbd2_cvb
