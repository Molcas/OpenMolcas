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

subroutine ddrestart2_cvb(c,axc,vec,hp,solp,maxdav,n,nvguess1,nvrestart1,eigval,eigvec)

use casvb_global, only: ifollow, nroot
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: maxdav, n, nvguess1, nvrestart1
real(kind=wp) :: c(n,maxdav), axc(n,maxdav), vec(n), hp(maxdav,maxdav), solp(maxdav), eigval(maxdav), eigvec(maxdav,maxdav)
integer(kind=iwp) :: ir, ir_use

call fmove_cvb(hp,eigvec,maxdav*maxdav)
call mxdiag_cvb(eigvec,eigval,maxdav)

nvrestart1 = 0
nvguess1 = 0
call mxatb_cvb(c,solp,n,maxdav,1,vec)
if (ifollow <= 2) then
  ! (Put lower-lying solutions in AxC :)
  do ir=1,nroot-1
    if (ifollow == 1) then
      ir_use = maxdav-ir+1
    else
      ir_use = ir
    end if
    call mxatb_cvb(c,eigvec(1,ir_use),n,maxdav,1,axc(1,ir+1))
  end do
  call fmove_cvb(axc(1,2),c(1,2),n*(nroot-1))
end if
call fmove_cvb(vec,c(1,1),n)

return

end subroutine ddrestart2_cvb
