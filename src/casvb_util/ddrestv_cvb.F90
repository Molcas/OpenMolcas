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

subroutine ddrestv_cvb(vec,avec,svec,ndim,ioffs,ause,suse)

use casvb_global, only: axc, c, maxd, nparm, nvguess, nvrestart, sxc
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ndim, ioffs
real(kind=wp) :: vec(ndim), avec(ndim), svec(ndim)
logical(kind=iwp) :: ause, suse

nvguess = nvguess+1
nvrestart = nvrestart+1
if ((nvguess > maxd) .or. (nvrestart > maxd)) then
  write(u6,*) ' Too many guess vectors in Davidson!',nvguess,nvrestart,maxd
  call abend_cvb()
end if
if (ndim+ioffs > nparm) then
  write(u6,*) ' Illegal call to DDRESTV :',ndim,ioffs,nparm
  call abend_cvb()
end if
call fzero(c(:,nvrestart),ioffs)
call fmove_cvb(vec,c(ioffs+1:,nvrestart),ndim)
call fzero(c(ndim+ioffs+1:,nvrestart),nparm-ioffs-ndim)
if (ause) then
  call fzero(axc(:,nvrestart),ioffs)
  call fmove_cvb(avec,axc(ioffs+1:,nvrestart),ndim)
  call fzero(axc(ndim+ioffs+1:,nvrestart),nparm-ioffs-ndim)
end if
if (suse) then
  call fzero(sxc(:,nvrestart),ioffs)
  call fmove_cvb(svec,sxc(ioffs+1:,nvrestart),ndim)
  call fzero(sxc(ndim+ioffs+1:,nvrestart),nparm-ioffs-ndim)
end if

return

end subroutine ddrestv_cvb
