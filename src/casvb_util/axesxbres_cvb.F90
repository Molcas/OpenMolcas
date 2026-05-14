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

subroutine axesxbres_cvb( &
#                        define _CALLING_
#                        include "ddres_interface.fh"
                        )

use casvb_global, only: ifollow, nroot
use Definitions, only: wp, iwp, u6

implicit none
#include "ddres_interface.fh"
real(kind=wp) :: alfa
integer(kind=iwp) :: i, nnegeig, nposeig

if (ifollow == 1) then
  nposeig = nroot-1
  nnegeig = itdav-nposeig
else if (ifollow == 2) then
  nnegeig = nroot-1
  nposeig = itdav-nnegeig
else
  nnegeig = 0
  nposeig = 0
  write(u6,*) ' Error in IFOLLOW with direct Fletcher!',ifollow
  call abend_cvb()
end if
res(:) = rhs(:)
do i=1,itdav
  if (i <= nnegeig) then
    alfa = eig_res
  else
    alfa = -eig_res
  end if
  res(:) = res(:)+(axc(:,i)-alfa*sxc(:,i))*solp_res(i)
end do

is_converged = .true.

return

end subroutine axesxbres_cvb
