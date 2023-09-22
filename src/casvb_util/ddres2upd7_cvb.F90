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

subroutine ddres2upd7_cvb(res,c,n)

use casvb_global, only: n_div
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: res(n), c(n)
real(kind=wp) :: resnrm1, resnrm2
real(kind=wp), external :: dnrm2_

if (n_div == 0) then
  call fmove_cvb(res,c,n)
else
  resnrm1 = dnrm2_(n_div-1,res(2),1)
  resnrm2 = dnrm2_(n-n_div,res(n_div+1),1)
  if (resnrm1 > resnrm2) then
    call fmove_cvb(res,c,n_div)
    call fzero(c(n_div+1),n-n_div)
  else
    call fzero(c,n_div)
    c(1) = res(1)
    call fmove_cvb(res(n_div+1),c(n_div+1),n-n_div)
  end if
end if
call ddproj_cvb(c,n)

return

end subroutine ddres2upd7_cvb
