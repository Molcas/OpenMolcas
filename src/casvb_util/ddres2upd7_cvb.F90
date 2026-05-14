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

subroutine ddres2upd7_cvb( &
#                         define _CALLING_
#                         include "ddres2upd_interface.fh"
                         )

use casvb_global, only: n_div
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "ddres2upd_interface.fh"
real(kind=wp) :: resnrm1, resnrm2
real(kind=wp), external :: dnrm2_

if (n_div == 0) then
  c(:) = res(:)
else
  resnrm1 = dnrm2_(n_div-1,res(2),1)
  resnrm2 = dnrm2_(n-n_div,res(n_div+1:),1)
  if (resnrm1 > resnrm2) then
    c(1:n_div) = res(1:n_div)
    c(n_div+1:) = Zero
  else
    c(1) = res(1)
    c(2:n_div) = Zero
    c(n_div+1:) = res(n_div+1:)
  end if
end if
call ddproj_cvb(c,n)

return

end subroutine ddres2upd7_cvb
