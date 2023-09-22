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

subroutine axbres_cvb(dum,sxc,rhs,res,solp_res,maxdav,n,itdav,eig_res,is_converged)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: maxdav, n, itdav
real(kind=wp) :: dum, sxc(n,maxdav), rhs(n), res(n), solp_res(maxdav), eig_res
logical(kind=iwp) :: is_converged

call mxatb_cvb(sxc,solp_res,n,itdav,1,res)
call dscal_(n,-eig_res,res,1)
call addvec(res,res,rhs,n)

is_converged = .true.

return
! Avoid unused argument warnings
if (.false.) call Unused_real(dum)

end subroutine axbres_cvb
