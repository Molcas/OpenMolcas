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

subroutine ddguess_cvb(vec,ndim,ioffs)

use casvb_global, only: c, maxd, nparm, nvguess
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndim, ioffs
real(kind=wp), intent(in) :: vec(ndim)

nvguess = nvguess+1
if (nvguess > maxd) then
  write(u6,*) ' Too many guess vectors in Davidson!',nvguess,maxd
  call abend_cvb()
end if
if (ndim+ioffs > nparm) then
  write(u6,*) ' Illegal call to DDGUESS :',ndim,ioffs,nparm
  call abend_cvb()
end if
c(1:ioffs,nvguess) = Zero
c(ioffs+1:ioffs+ndim,nvguess) = vec(:)
c(ioffs+ndim+1:,nvguess) = Zero

return

end subroutine ddguess_cvb
