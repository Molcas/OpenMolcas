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

subroutine o7b_cvb( &
#                  define _CALLING_
#                  include "optb_interface.fh"
                  )

use casvb_global, only: expct, have_solved_it, hh, ip, odx, scalesmall
use casvb_interfaces, only: ddasonc_sub, ddres2upd_sub
use Constants, only: One, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "optb_interface.fh"
integer(kind=iwp) :: ioptc2, ipu, iter2
real(kind=wp) :: fac1, fx_exp, resthr_old = -One, resthr_use
logical(kind=iwp) :: skip
real(kind=wp), external :: dnrm2_
procedure(ddasonc_sub) :: asonc7_cvb
procedure(ddres2upd_sub) :: ddres2upd7_cvb

if (.not. close2conv) then
  resthr_use = 1.0e-5_wp
else
  resthr_use = 0.05_wp*grdnrm
  resthr_use = min(1.0e-5_wp,resthr_use)
  resthr_use = max(1.0e-9_wp,resthr_use)
end if
skip = ((resthr_use == resthr_old) .and. have_solved_it)
resthr_old = resthr_use
if (.not. skip) then
  call axex_cvb(asonc7_cvb,ddres2upd7_cvb,odx,resthr_use,ioptc2,iter2,fx_exp)
  have_solved_it = .true.
  expct = Half*fx_exp

  if (ip >= 2) write(u6,'(a,i4)') ' Number of iterations for direct diagonalization :',iter2

  if (ioptc2 /= 0) then
    write(u6,*) ' Direct diagonalization not converged!'
    call abend_cvb()
  end if

  if (ip >= 2) then
    write(u6,'(a)') ' Eigenvector to be followed :'
    call vecprint_cvb(odx,nparm+1)
  end if
  if (abs(odx(1)) > 1.0e-8_wp) then
    fac1 = one/odx(1)
  else
    fac1 = sign(one,odx(1))
  end if
  odx(1:nparm) = fac1*odx(1:nparm)
  odx(1:nparm) = odx(2:nparm+1)
end if
dxnrm = dnrm2_(nparm,odx,1)
if (.not. close2conv) then
  ipu = 1
else
  ipu = 2
end if
if ((dxnrm > hh) .or. scalesmall(ipu)) then
  odx(1:nparm) = hh/dxnrm*odx(1:nparm)
  dxnrm = hh
end if

return

end subroutine o7b_cvb
