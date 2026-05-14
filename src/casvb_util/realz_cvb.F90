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

subroutine realz_cvb(arr,nmax,nread,ifc)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nmax, ifc
real(kind=wp), intent(out) :: arr(nmax)
integer(kind=iwp), intent(out) :: nread
integer(kind=iwp), parameter :: nbuf = 100
integer(kind=iwp) :: nleft, nread1
real(kind=wp) :: tmp(nbuf)

nread = 0
do
  tmp(:) = Zero
  nleft = nmax-nread
  call real_cvb(tmp,min(nbuf,nleft),nread1,ifc)
  arr(nread+1:nread+nread1) = tmp(1:nread1)
  nread = nread+nread1
  if (nread1 <= 0) exit
end do

return

end subroutine realz_cvb
