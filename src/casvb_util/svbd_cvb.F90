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

subroutine svbd_cvb(orbs,cvb,fx,ioptc,iter)

use casvb_global, only: civb1, civb2, cvbdet, gjorb2, ifollow, ipdd, ipr, isaddle, isaddledd, follow, have_solved_it, maxdav, &
                        memplenty, mxiter, n_div, norb, nortiter, nroot, nvb, orththr, ovraa, resthr, svb
use casvb_interfaces, only: ddasonc_sub, ddres_sub, ddres2upd_sub, ddsol_sub
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: orbs(norb,norb)
real(kind=wp), intent(inout) :: cvb(nvb)
real(kind=wp), intent(out) :: fx
integer(kind=iwp), intent(out) :: ioptc, iter
integer(kind=iwp) :: ifollow1, nvguess, nvrestart
real(kind=wp), allocatable :: c(:,:), dum(:,:), res(:), rhs(:), rhsp(:), solp(:), solp_res(:), sxc(:,:)
procedure(ddasonc_sub) :: asonc1_cvb
procedure(ddres_sub) :: ddressvb_cvb
procedure(ddres2upd_sub) :: ddres2upd10_cvb
procedure(ddsol_sub) :: ddsolsvb_cvb

call makegjorbs_cvb(orbs)

if (memplenty) then
  call cicopy_cvb(civb1,civb2)
else
  call cird_cvb(civb2,61001.2_wp)
end if
call applyt_cvb(civb2,gjorb2)
call ci2vbg_cvb(civb2,cvbdet)
call mma_allocate(rhs,nvb,label='rhs')
call vb2strg_cvb(cvbdet,rhs)

ioptc = 1
nvguess = 1
nvrestart = 0
call mma_allocate(c,nvb,maxdav,label='c')
c(:,1) = cvb(:)
if (.not. follow) then
  ifollow1 = 1
else
  ifollow1 = 4
end if

ifollow = ifollow1
isaddledd = isaddle
nroot = max(1,isaddledd+1)
ipdd = ipr(3)
n_div = 0
call mma_allocate(sxc,nvb,maxdav,label='sxc')
call mma_allocate(res,nvb,label='res')
call mma_allocate(rhsp,maxdav,label='rhsp')
call mma_allocate(solp,maxdav,label='solp')
call mma_allocate(solp_res,maxdav,label='solp_res')
call mma_allocate(dum,max(nvb,maxdav),maxdav,label='dum')
call dirdiag_cvb(asonc1_cvb,ddsolsvb_cvb,ddressvb_cvb,ddres2upd10_cvb,c,dum,sxc,.false.,cvb,res,rhs,dum,rhsp,solp,solp_res, &
                 .false.,.false.,.true.,maxdav,nvb,nvb,nvguess,nvrestart,isaddle,ifollow1,mxiter,resthr,orththr,nortiter,Zero, &
                 ioptc,iter,fx,ipr(3))
call mma_deallocate(c)
call mma_deallocate(sxc)
call mma_deallocate(res)
call mma_deallocate(rhs)
call mma_deallocate(rhsp)
call mma_deallocate(solp)
call mma_deallocate(solp_res)
call mma_deallocate(dum)
have_solved_it = .true.

ovraa = One
svb = fx

return

end subroutine svbd_cvb
