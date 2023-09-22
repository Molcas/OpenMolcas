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

subroutine svbd2_cvb(orbs,cvb,fx,ioptc,iter,civec,civbs,gjorb,gjorb2,gjorb3,cvbdet,c,sxc,res,rhs,rhsp,solp,solp_res)

use casvb_global, only: follow, have_solved_it, nortiter, orththr, resthr
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
#include "optze_cvb.fh"
real(kind=wp) :: orbs(norb,norb), cvb(nvb), fx, civec(ndet), civbs(ndet), gjorb(*), gjorb2(*), gjorb3(*), cvbdet(ndetvb), &
                 c(nvb,maxdav), sxc(nvb,maxdav), res(nvb), rhs(nvb), rhsp(maxdav), solp(maxdav), solp_res(maxdav)
integer(kind=iwp) :: ioptc, iter
#include "print_cvb.fh"
integer(kind=iwp) :: ifollow, nvguess, nvrestart
real(kind=wp) :: dum(max(nvb,maxdav),maxdav) !IFG
external :: asonc1_cvb, ddres2upd10_cvb, ddressvb_cvb, ddrestart_cvb, ddsolsvb_cvb

call makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)

if (memplenty) then
  call cicopy_cvb(civec,civbs)
else
  call cird_cvb(civbs,61001.2_wp)
end if
call applyt_cvb(civbs,gjorb2)
call ci2vbg_cvb(civbs,cvbdet)
call vb2strg_cvb(cvbdet,rhs)

ioptc = 1
nvguess = 1
nvrestart = 0
call fmove_cvb(cvb,c,nvb)
if (.not. follow) then
  ifollow = 1
else
  ifollow = 4
end if

call ddinitsvb_cvb(ifollow,isaddle,ip(3))
call ddres2updinit_cvb(0)
call dirdiag_cvb(asonc1_cvb,ddsolsvb_cvb,ddressvb_cvb,ddres2upd10_cvb,ddrestart_cvb,c,dum,sxc,.false.,cvb,res,rhs,dum,rhsp,solp, &
                 solp_res,.false.,.false.,.true.,maxdav,nvb,nvb,nvguess,nvrestart,isaddle,ifollow,mxiter,resthr,orththr,nortiter, &
                 Zero,ioptc,iter,fx,ip(3))
have_solved_it = .true.

ovraa = one
svb = fx

return

end subroutine svbd2_cvb
