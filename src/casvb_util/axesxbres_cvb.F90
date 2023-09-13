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

subroutine axesxbres_cvb(axc,sxc,rhs,res,solp_res,maxdav,n,itdav,eig_res,is_converged)

implicit real*8(a-h,o-z)
logical is_converged
#include "direct_cvb.fh"
dimension axc(n,maxdav), sxc(n,maxdav), rhs(n), res(n)
dimension solp_res(maxdav)

if (ifollow == 1) then
  nposeig = nroot-1
  nnegeig = itdav-nposeig
else if (ifollow == 2) then
  nnegeig = nroot-1
  nposeig = itdav-nnegeig
else
  nnegeig = 0
  nposeig = 0
  write(6,*) ' Error in IFOLLOW with direct Fletcher!',ifollow
  call abend_cvb()
end if
call fmove_cvb(rhs,res,n)
do i=1,itdav
  if (i <= nnegeig) then
    alfa = eig_res
  else
    alfa = -eig_res
  end if
  do ivb=1,n
    res(ivb) = res(ivb)+(axc(ivb,i)-alfa*sxc(ivb,i))*solp_res(i)
  end do
end do

is_converged = .true.

return

end subroutine axesxbres_cvb
