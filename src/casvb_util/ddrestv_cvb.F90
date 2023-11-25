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
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndim, ioffs
real(kind=wp), intent(in) :: vec(ndim), avec(ndim), svec(ndim)
logical(kind=iwp), intent(in) :: ause, suse

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
c(1:ioffs,nvrestart) = Zero
c(ioffs+1:ioffs+ndim,nvrestart) = vec(:)
c(ioffs+ndim+1:,nvrestart) = Zero
if (ause) then
  axc(1:ioffs,nvrestart) = Zero
  axc(ioffs+1:ioffs+ndim,nvrestart) = avec(:)
  axc(ioffs+ndim+1:,nvrestart) = Zero
end if
if (suse) then
  sxc(1:ioffs,nvrestart) = Zero
  sxc(ioffs+1:ioffs+ndim,nvrestart) = svec(:)
  sxc(ioffs+ndim+1:,nvrestart) = Zero
end if

return

end subroutine ddrestv_cvb
