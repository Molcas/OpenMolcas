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

subroutine o5b2_cvb(nparm,dx,grad,dxnrm,close2conv)

use casvb_global, only: hh, maxize, scalesmall
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nparm
real(kind=wp) :: dx(nparm), grad(nparm), dxnrm
logical(kind=iwp) :: close2conv
integer(kind=iwp) :: ipu
real(kind=wp), external :: dnrm2_

call fmove_cvb(grad,dx,nparm)
if (.not. maxize) call dscal_(nparm,-One,dx,1)
dxnrm = dnrm2_(nparm,dx,1)
if (.not. close2conv) then
  ipu = 1
else
  ipu = 2
end if
if ((dxnrm > hh) .or. scalesmall(ipu)) then
  call dscal_(nparm,hh/dxnrm,dx,1)
  dxnrm = hh
end if

return

end subroutine o5b2_cvb
