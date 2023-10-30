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

subroutine axex_cvb(asonc,ddres2upd,vec,resthr_inp,ioptc,iter,fx_exp)

use casvb_global, only: c, nparm, solp
use casvb_interfaces, only: ddasonc_sub, ddres2upd_sub
use Definitions, only: wp, iwp

implicit none
procedure(ddasonc_sub) :: asonc
procedure(ddres2upd_sub) :: ddres2upd
real(kind=wp), intent(out) :: vec(nparm), fx_exp
real(kind=wp), intent(in) :: resthr_inp
integer(kind=iwp), intent(out) :: ioptc, iter

call axesx2_cvb(asonc,ddres2upd,vec,resthr_inp,ioptc,iter,fx_exp,c,.true.,solp)

return

end subroutine axex_cvb
