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

subroutine ddrhs_cvb(vec,ndim,ioffs)

use casvb_global, only: mxrhs, nparm, nvrhs, rhs
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndim, ioffs
real(kind=wp), intent(in) :: vec(ndim)

nvrhs = nvrhs+1
if (nvrhs > mxrhs) then
  write(u6,*) ' Too many RHS vectors in Davidson!',nvrhs,mxrhs
  call abend_cvb()
end if
if (ndim+ioffs > nparm) then
  write(u6,*) ' Illegal call to DDRHS :',ndim,ioffs,nparm
  call abend_cvb()
end if
rhs(1:ioffs,nvrhs) = Zero
rhs(ioffs+1:ioffs+ndim,nvrhs) = vec(:)
rhs(ioffs+ndim+1:,nvrhs) = Zero

return

end subroutine ddrhs_cvb
