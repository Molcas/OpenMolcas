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

subroutine ddres7_cvb(axc,sxc,dum,res,solp_res,maxdav,n,itdav,eig_res,is_converged)

use casvb_global, only: iroot, jroot
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: maxdav, n, itdav
real(kind=wp) :: axc(n,maxdav), sxc(n,maxdav), dum, res(n), solp_res(maxdav), eig_res
logical(kind=iwp) :: is_converged
integer(kind=iwp) :: i, ivb

call fzero(res,n)
do i=1,itdav
  do ivb=1,n
    res(ivb) = res(ivb)+(axc(ivb,i)-eig_res*sxc(ivb,i))*solp_res(i)
  end do
end do

is_converged = (jroot == iroot)

return
! Avoid unused argument warnings
if (.false.) call Unused_real(dum)

end subroutine ddres7_cvb
