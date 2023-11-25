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

subroutine o10b_cvb( &
#                   define _CALLING_
#                   include "optb_interface.fh"
                   )

use casvb_global, only: have_solved_it, ip, odx
use casvb_interfaces, only: ddasonc_sub, ddres2upd_sub
use Definitions, only: wp, iwp, u6

implicit none
#include "optb_interface.fh"
integer(kind=iwp) :: ioptc2, iter2
real(kind=wp) :: fx_exp, resthr_use
real(kind=wp), external :: dnrm2_
procedure(ddasonc_sub) :: asonc10_cvb
procedure(ddres2upd_sub) :: ddres2upd10_cvb

if (.not. close2conv) then
  resthr_use = 1.0e-5_wp
else
  resthr_use = 0.05_wp*grdnrm
  resthr_use = min(1.0e-5_wp,resthr_use)
  resthr_use = max(1.0e-9_wp,resthr_use)
end if
call axexb_cvb(asonc10_cvb,ddres2upd10_cvb,odx,resthr_use,ioptc2,iter2,fx_exp)
have_solved_it = .true.

if (ip >= 2) write(u6,'(a,i4)') ' Number of iterations for direct diagonalization :',iter2

if (ioptc2 /= 0) then
  write(u6,*) ' Direct diagonalization not converged!'
  call abend_cvb()
end if

dxnrm = dnrm2_(nparm,odx,1)

return

end subroutine o10b_cvb
